# Biological interpretation claim audit (#179)

Audit date: 2026-09-14 (updated 2026-10-07). Policy: [shared interpretation guide](biological_interpretation.md).
This is a wording/provenance audit, not a reanalysis or validation of mechanisms.
`Corrected here` means the source prose/template on this branch is changed;
previously exported images still require regeneration. `Proposed correction`
records reviewable wording for another branch, not an assertion that its owner
has applied it. `Retain` means the scoped wording/evidence limit is appropriate.

## Inspected scope and reproducibility

| Surface | Inspected revision and coverage | Responsible adoption |
| --- | --- | --- |
| Current checkout | `c9b51a7` / `3b73376`: README, CONTRIBUTING, AGENTS, generated Copilot instructions, PR template, all methods, docs, profile aggregation and enrichment/depletion scripts | All contributors; Noel for current guild/Table 1 methods and figures; manuscript/public-description owners to be identified by coordinators |
| Proposed overlay, PR #173 | `85f97ac714c6d7c580c3a548b7f428229d58b634`: Table 1 methods, `guild_maps.py::_significance_note`, overlay source and output schema; figure captions inspected in generator source | Noel, with Jed and Ilyas for their statistical labels |
| Rosalind archive/review, PR #192 | `593f3c88100c915187c8bc06e2cbdbedef7730c3`: `docs/rosalind-six-trait-interpretation.md`; also supplied in the issue request. The archived review describes its own evidence at `f22b31bd4bdef6e6e55bb5830ca631b81e3d4319` | All narrative authors; follow-up leads below |
| Taxon contribution/influence audit, PR #194 | `021e5b449080c4a4f2f22a2ae329fccffc7eb367`: `methods/table1_taxon_influence.md` and `methods/table1_taxon_influence_results.md`; inspected the measured contribution, concentration, shared-fraction and bounded lineage-exclusion summaries | Jed with Noel and Ilyas; retain the provisional annotation-level limits and integrate the missing native results before closing #190 |
| Composition branch | `a2b2df98e9beec98b772c9cb5af7aa090d033a98`: targeted review of generated narrative/title sites in `analyze_parkinsons_disease_trait_profiles.py`, `investigate_facultative_aerobe_expansion.py`, `investigate_gut_anaerobe_depletion.py`, and figure annotation templates | Jed; remaining branch artifact owners to inventory |
| Enrichment branch | `c9b51a7`: `scripts/generate_enrichment_summary.py`, `scripts/generate_depletion_summary.py`, `docs/milestone_taxon_model_before_after_summary.md`, and `docs/draft_pr_overview_and_merge_guide.md` | Ilyas |

The PR revisions were fetched through `refs/pull/173/head` and
`refs/pull/192/head`; no changes were made to those PRs. Other branches were
inspected with `git show`/`git grep`. Reproduce a claim lookup with
`git show <revision>:<path>` and the quoted phrase below. The separate-branch
reviews are targeted, not exhaustive certification of every historical image.
Unpublished manuscripts, slides and public pages not present in this checkout
remain to be inventoried by Jed, Ilyas and Noel, with a named owner and revision.
The source evidence audits and validation work in #188-#191 remain open work.

## Current prose, workflow and generated captions

| ID / location | Affected claim or omission | Correction / evidence and disposition |
| --- | --- | --- |
| C01 `methods/guild_condition_body_site_heatmap.md`, Cell states and row selection | "no taxon ... carried" the guild; "guild is simply absent" | **Corrected here:** no recorded guild annotation signal under selected sources/mappings. Distinguish unsampled, low support and zero signal; none establishes biological absence. |
| C02 `src/bsdbng/guild_maps.py`, `write_figure` and `write_body_site_panels` captions/docstrings | Dark grey means "no taxon carried" the trait; "an absent trait" | **Corrected here:** no recorded annotation signal; captions retain annotation-share units and state the evidence limit and curated meaning of direction. Split long caption/coverage lines so these limits remain visible in the rendered panels. |
| C03 `methods/table1_trait_condition_body_site_maps.md`, cell states, ties and empty sides | "trait is absent from both sides"; profiles whose taxa "carry none" | **Corrected here:** zero recorded counts, not biological absence; empty-side handling is an analytical policy. |
| C04 same methods, right-panel correlations; `table1_trait_maps.py` comments/docstring | Contrast correlations "do not share that artefact"; coccus/sphere are "close to the same measurement" | **Corrected here:** shared denominators, taxa and selection can induce dependence even across paired contrasts. Correlation does not establish semantic equivalence, independent mechanisms or justify merging IDs. #165/#182 and #189 retain their work. |
| C05 `methods/three_arm_statistics.md`, Compositionality and Vintage | Any pair of fractions is negatively correlated "by construction"; paired tests "do not share the artefact"; same dump implies valid comparison | **Corrected here:** closure constrains dependence for positive denominators but does not force every pair negative; paired contrasts remain denominator-sensitive. Match populations, denominator, transformation, background and family, not just vintage. |
| C06 `table1_trait_maps.py::write_shift_figure`, title, subtitle and axis | "directional shift" without an explicit annotation qualifier on the standalone figure | **Corrected here:** paired annotation contrast; show study support, ties, intervals, q-family and curated direction, with no patient measurement claim. Numerical computation is unchanged. |
| C07 README, methods index and per-figure methods | Public project scope and per-figure descriptions lack a common interpretation contract | **Corrected here:** link the guide/audit and state annotation scope, contrast-orientation restriction and descriptive heatmap uncertainty limits. |
| C08 CONTRIBUTING, AGENTS, Copilot instructions, PR checklist | Review workflow does not require biological claim/provenance review | **Corrected here:** guide applies to every author/reviewer and output; require effect/uncertainty/support/inference/scope, caption/template review and audit updates. Regenerate Copilot via `just sync-agent-docs`. |
| C09 Table 1 methods, agreement with composition arm | Explicitly describes shared-matrix agreement as robustness rather than independent evidence | **Retain**, with denominator/selection/calibration qualifications from the guide. Do not promote agreement to replication. |
| C10 `MANUSCRIPT_AUDIT.md`, #131/#132/#133 and uncertainty notes | Existing warnings about assertion-level inference and incomplete manuscript access | **Retain:** where tracked, link this biological audit. A statistical audit or corrected implementation does not validate biological interpretations. |
| C11 `docs/origin-story.md`, "taxa and their abundance directions"; access/benchmark docs | Ingestion history describes the original curated direction; access docs describe data acquisition | **Retain** as provenance descriptions, not claims of patient-level trait measurement. No mechanistic result claim identified in those documents. |

## PR #173: proposed overlay and narrative

Source: [PR #173](https://github.com/CultureBotAI/microbial-trait-profiling/pull/173),
revision above, `methods/table1_trait_condition_body_site_maps.md`, "Third layer"
and "What the outlines show"; `src/bsdbng/guild_maps.py::_significance_note`.
All entries below are **proposed corrections for the PR owner**, not edits to
the separate branch. Reported numbers are historical values in that PR and have
not been recomputed or approved here.

| ID / affected claim | Replacement or supporting evidence | Owner / remaining work |
| --- | --- | --- |
| P01 Both arms answer whether a cell's share is "more than chance" | The outline summarizes profile-level annotation enrichment calls against a specified corpus background. A triangle reports a body-site pooled paired composition test. Neither supplies a common cell-level test. | Noel with Ilyas/Jed; state each null, unit, family, background and display aggregation rule. |
| P02 Black outline: studies "are significantly enriched ... on their own" | For enrichment-derived marks: "Outline: a strict majority of studies contributing to this cell have qualifying profile-level annotation calls against the filtered corpus (at least two total studies; enrichment uses profile q <= 0.001 and fold >= 1.5). This aggregation has no separately calibrated cell-level p/q." State which verdicts actually contribute; see P10 for the implementation/display mismatch. The figure's masking still applies. | Noel/Ilyas; attach called/eligible counts, observed/expected/fold values and profile-to-study aggregation; #131/#182 limitations remain. |
| P03 Triangle points to the side "carrying more of it" | "Triangle: significant increased-minus-decreased transformed annotation-share contrast pooled across this body site (q <= 0.05 in its stated family); points to the higher-share direction. It does not test each condition in the column." | Noel/Jed; carry transformed units, effect, uncertainty and study counts. |
| P04 Fecal six-trait pattern is "the oxygen hypothesis stated in one line"; organisms rising in "a dysbiotic gut" are oxygen-tolerant | "In the retained fecal corpus, the reported transformed annotation contrasts favor the increased signatures for five traits and decreased signatures for cellulose catabolic process. These associations motivate an intestinal redox/substrate hypothesis; they do not establish oxygen exposure, functional activity or disease causation." | Noel; verify disease/control orientation (#188), trait semantics (#189) and actual measurements (#191). Do not translate all five traits into oxygen tolerance. |
| P05 Only one contrast survives outside the gut | "Under these filters/tests, only one non-fecal contrast passes the stated threshold; other effects and their uncertainty/support must also be reported." A non-significant or unavailable result is not absence of a biological effect. | Noel; retain #118/#120/#176 stratum scope and #183 sensitivity. |
| P06 Chronic constipation outline is the "opposite of the fecal contrast" and what "slow transit favouring fibre degraders would look like" | "Cellulose-associated annotations meet the enrichment display rule in increased constipation signatures. Enrichment against a background is not the same estimand as the pooled paired fecal contrast. Slow transit is a hypothesis requiring measured transit/substrate data and a within-stratum contrast." | Noel with #188/#191; an outline alone does not prove a reversed paired effect. |
| P07 Cervical "direction inverts"; "vaginal community is rod-dominated at baseline" | "Coccus-associated annotations meet the enrichment display rule in decreased signatures for the named cervical/vaginal strata. This is not a test of a reversed paired effect or of baseline community morphology." A baseline claim needs a scoped source and matched observations. | Noel; #188/#189, exact anatomical/condition scope; do not generalize a gut mechanism. |
| P08 Named condition patterns are a "coherent set"; pattern "crosses" feces and saliva | Describe the specific annotation calls and per-cell study support. Shared studies, traits and corpus selection do not make these independent confirmation or multiple independent mechanisms. | All three coordinators; carry overlapping-study IDs and #165/#182 limitations. |
| P09 Outline means "unusually assertion-rich", not "taxon-rich"; depends on unresolved urn issue | **Retain** that limit, but explicitly tie "unusually" to the specified background/test and calibration. A fixes-in-progress warning is not evidence that the biological explanation is correct. | Ilyas/Noel; #131/#132/#182 remain independent requirements. |
| P10 `trait_significance.py::build_overlay` versus `_significance_note` | The code counts studies with any nonempty verdict toward the majority, labels enrichment/depletion by their relative counts, and can mark a cell for condition-level composition alone or both arms. The legend describes all outlines as enrichment. Require a legend keyed to the actual mark kind and verdict, or restrict marks to the declared rule; report any depletion threshold and its test. Do not label every outline as an enrichment majority. | Noel with Ilyas/Jed; reconcile the generated legend and actual mark provenance in #172/#173 before interpreting the 62-cell summary. |

## PR #192: worked claim corrections from the separate review

Source: [archived transcript and repository review](https://github.com/CultureBotAI/microbial-trait-profiling/pull/192).
The transcript is intentionally preserved as supplied. **Retain the archive and
its review; use the corrections below in new narratives.** No new biological
evidence is supplied by this issue.

| ID / archived claim | Required correction/evidence | Follow-up |
| --- | --- | --- |
| R01 "probably causal ecological readouts" | Unsupported. Annotation contrasts motivate falsifiable hypotheses about redox, substrate availability and taxonomic turnover; require independent measurements and a causal design for causal language. | #179 guide; #191 validation feasibility |
| R02 Disease pushes all host-associated microbiomes toward one gut-redox state | Limit the candidate mechanism to selected intestinal settings; report exact anatomy, eligible contrast, effects, uncertainty and retained support. Human/paired filtering does not establish disease versus healthy control. | #188 Noel with Jed; #118/#120/#176 |
| R03 Veillonella is "facultative-adjacent" evidence for oxygen use | The review identifies anaerobic nitrate respiration in Rojas-Tapias et al. (2022), DOI 10.1038/s41564-022-01224-7. This is not evidence of facultative aerobiosis for this annotation. | #189 Ilyas with Jed |
| R04 Coccus/sphere should be collapsed; effect stays the same | Distinct IDs METPO:1000668 / METPO:1000683; distinguish ontology relationships, mapping identity and empirical overlap. Require source assertions and a deduplicated grouping with separately defined denominator/testing family. Preserve original results. | #189; calibration #182 |
| R05 Arsenate feature means arsC, plasmids and co-selection | The qualified pathway feature does not identify a gene or mobile element. Classify source evidence as detoxification/respiratory/mixed/unresolved; do not infer arsC dominance or arrAB absence. | #189 |
| R06 Moderately halophilic means halotolerant/osmotolerant; disease creates salt stress | Recover source growth requirement/optimum/range, units and culture conditions for METPO:1000623; measurements are needed for host salt exposure. Do not infer genera from the label. | #189; measured evidence #191 |
| R07 Cellulose catabolism implies lost butyrate production/cross-feeding | GO:0030245 is not measured complete degradation, SCFA synthesis, expression or flux. Distinguish those evidence levels. Same-annotation correlations are not independent validation. | #191 Noel with Ilyas; #165/#182 |
| R08 Named opportunists should dominate; exclusions should weaken effects | Reconcile all six traits and both directions with source/rollup counts. Predeclare the taxonomic universe and exclusions; preserve opposite/null findings and recompute denominators/backgrounds/eligibility. An exclusion is an influence diagnostic. | #190 Jed with Noel/Ilyas |
| R09 Inflammation and SCFA correlations will validate the model | Inventory actual measured markers at the usable observation level, covariates and cohort overlap; define a feasible independent protocol before testing. Disease labels and exclusion flags are not measured inflammation or exposure. | #191 |

The PR review's cited primary studies support mechanisms in their own settings;
they do not explain these six repository directions by themselves. Its finding
of missing upstream KG artifacts applies to its inspected revision/environment,
not a claim that the original publications lack richer evidence.

## Other generated summaries, scripts, and documentation

| ID / source | Affected claim | Replacement and required evidence | Owner / Status |
| --- | --- | --- | --- |
| G01 `scripts/generate_enrichment_summary.py`, `summary_text` | "clinical functional burden independent of cohort size" | "Mean number of annotation-enrichment calls per eligible profile under the specified background and thresholds." Neither clinical burden nor independence from cohort size follows from profile normalization. | Ilyas; **Corrected here** |
| G02 `scripts/generate_enrichment_summary.py`, Panel C/text headers | "Anatomical Enrichment", "Hits/100 Samples" | "Anatomical Annotation Representation", "Calls / 100 Reported Samples". Normalized calls relative to reported sample totals, with note that this is an aggregate call ratio across signatures, not sample-level trait frequency. | Ilyas; **Corrected here** |
| G03 `scripts/generate_depletion_summary.py`, `summary_text` | "Functional Blind Spots", "average functional loss", traits "avoided in disease signatures", "functional depletion" | "Under-Represented Trait Scope", "mean under-representation calls per profile against the specified background (zero mapped counts do not establish biological loss or absence)." | Ilyas; **Corrected here** |
| G04 `scripts/generate_depletion_summary.py`, Panel C/colorbar and text output | "Anatomical Depletion", "Losses / 100 Samples" | "Anatomical Under-Representation", "Calls / 100 Reported Samples". Aggregate call ratio across signatures, not patient-level functional losses. | Ilyas; **Corrected here** |
| G05 `analyze_parkinsons_disease_trait_profiles.py`, title and generated narrative | "PD feces anaerobe-depletion target taxa"; "Anaerobe depletion profile-level effect" | "Anaerobe-associated annotation contrasts in retained Parkinson's disease-labeled fecal signatures," or "taxon signature occurrence counts" where that is the plotted quantity. Preserve its existing non-causal caveat and supply orientation/support before disease-level wording. | Jed; Proposed correction (#188/#190) |
| G06 expansion/depletion investigation script titles, output slugs and figure annotation templates | "aerobe expansion" / "anaerobe depletion" can read as biological measurements | Preserve historical filenames/IDs for traceability; qualify display titles and generated narratives with the actual annotation or signature-occurrence metric and source-defined direction. Candidate-taxon predictions remain hypotheses. | Jed; Proposed correction |
| G07 `trait-visualization/trait-profile-investigation/paired_contrast_interpretation_notes.md` | Notes distinguish z-scores from raw abundance and curated direction, and warn against combining dependent p-values | **Retain** these useful starting rules. Name the actual denominator/version, support and uncertainty; descriptive screening language does not remove calibration requirements. | Jed |
| G08 same notes, interpretation/integration grids | Positive mean z means "enriched"; single-trait enrichment treated as prevalence; "High enrichment" implies "Common trait" | Replace with "above-average transformed annotation signal" for z-scores. Define ORA as a comparison against its background; a large fold need not mean common/high-prevalence. Keep prevalence, background enrichment and paired contrast in separate columns. | Jed with Ilyas; Proposed correction |
| G09 `trait-visualization/trait-profile-investigation/six-shared-traits/interpretations.md`, Overall/Cross-Disease Patterns | Traits "are not disease-specific"; "shared microbial phenotype axis" including "aerotolerance/oxygen exposure" and "halophily/osmotic tolerance" | "Selected annotation signals recur across the displayed condition/site groups." Display recurrence does not test disease specificity, measure oxygen exposure or justify salt-label translation. Keep native semantics, eligible orientation and anatomical limits (#188/#189). | Jed; Proposed correction |
| G10 `trait-visualization/figure_methods.md`, taxon heatmap legend | Zero "cells indicate absence/no annotation" | "Zero cells indicate no recorded annotation under the selected mappings." Distinguish missing mapping/coverage and biological absence; update the generator and rendered legend too. | Jed; Proposed correction |
| G11 PR #194 contribution/influence results | Largest contributing genera, pooled baseline contrasts and largest changes in a predeclared 22-exclusion panel | **Retain as provisional annotation influence diagnostics:** carry resolved-rank coverage, exact corpus/denominator, pairs/studies, effect units, conditional intervals, ambiguity and selection limits. Count concentration is not an additive fraction-effect contribution, and an exclusion is not a causal intervention. Do not infer patient abundance, oxygen exposure, SCFA production or a mechanism. Orientation (#188), semantics/source identity (#189), native Jed/ORA integration and calibration (#182) remain unresolved. | Jed with Noel/Ilyas; PR #194 / #190 |
| G12 `docs/milestone_taxon_model_before_after_summary.md` | Claimed restored exchangeability, calibrated probabilities, genuine biological drivers, demonstrated recovery power, and specific disease/control trait rankings while referencing historical artifacts with p=0 rows | **Corrected here:** Replace unsupported completion and biological claims with measured implementation status. Separate intended taxon urn architecture from demonstrated empirical benchmark results. Clarify that binary presence alone does not establish independent sampling or universal biological exchangeability across phylogenetically correlated taxa. State that exact discrete hypergeometric p-values need not be uniform under the null. Describe recurrent annotation signals as representation in curated signatures rather than genuine biological drivers or disease/control rankings. Note that historical `c75744a` artifacts contained p=0 rows and empirical calibration benchmarking awaits #263 / #182. | Ilyas; **Corrected here** |
| G13 `docs/draft_pr_overview_and_merge_guide.md` | Claimed completion of calibrated probabilities, unified production adapter, and demonstrated recovery power across branch deliverables | **Corrected here:** Separate implemented infrastructure/interfaces from demonstrated calibration results. Note that calibration benchmarking requires execution of the full synthetic and empirical battery (#263 / #182) and that biological wording is governed by `methods/biological_interpretation.md`. | Ilyas; **Corrected here** |

## Adoption record

The repository guide, claim audit, workflow integration, generator corrections (G01-G04),
and milestone documentation updates (G12-G13) are delivered on #179. PR #173, PR #194's
pending native integrations, separate analysis branches, unpublished manuscript sections,
and public outputs still require their owners to apply the proposed corrections and
record artifact revisions. No coordinator sign-off or outside-branch modification is
asserted here. Add each newly identified output to this ledger; record its owner,
source revision, exact claim, evidence or replacement, and completion revision.
#188-#191 complement this guide and do not close biological validation.
