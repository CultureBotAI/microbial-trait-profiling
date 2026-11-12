# Presence and annotation-depth sensitivity (#772)

This additive analysis is declared before the coordinated rerun. It consumes a
verified core bundle, preserves its contrast population and already applied
direction orientation, and leaves the three arms' primary calls and plans
unchanged. It is an annotation sensitivity analysis, not a causal adjustment or
measurement of physiology.

Protocol `presence-depth-study-v2` follows review #892. The complete candidate
family is the core catalog's non-MediaDive vocabulary, in the actual profile
header order, including candidates with zero observed signal. Its actual size and
fingerprint are recorded in an `AnalysisBasis`; no historical feature count is
hard-coded. Each representation and test has a separately declared full family,
BY as the primary correction and BH as a comparison. Both the full screen and
the focal table retain these primary full-family values in `q_by` and `q_bh`,
with `family_size` identifying the full number of candidates. The focal table's
generic `family_id` and `hypothesis_id` name that full family and its exact
retained member. `full_family_id` and `full_hypothesis_id` are explicit aliases.
The additional six-trait family has separate `family_id_six_trait` and
`hypothesis_id_six_trait` identities; its results appear only in `q_by_six_trait`,
`q_bh_six_trait`, and `family_size_six_trait`. The six-trait BH comparison uses
the focal Table 1 membership and procedure within this sensitivity's declared
representation and test; it is not the same estimand or population as every
Table 1 result. Missing focal features retain unavailable slots.
Family declarations and correction records use the authoritative LinkML
multiple-testing v1.1 models and `mtp.multiple_testing`; they are replay-verified
against their declared inputs before writing. The shared policy also reports BY
for the six-trait family, without replacing the full-screen primary call.
Candidates with all ties or zero residual variance are untestable, as are
unsupported or insufficient fits. Untestable candidates occupy their family
slots but receive no reported p-value or q-value; internal
correction uses p=1 for these slots. Comparisons across these sensitivity families
are exploratory and are not independent replication.

## Focal export contract and historical migration

New `focal_traits.tsv` rows carry
`export_contract_version=presence-depth-focal-export-v2`. The same value is
recorded as `focal_export_contract_version` in `provenance.json` and stated in
the generated README. This export version is separate from scientific protocol
`presence-depth-study-v2`: no family definition, member inventory, p-value,
correction procedure, q-value or family size changes. The producer projects
copies of the focal rows only after the full and six-trait corrections have
completed; it does not rerun either correction or modify its input rows.

| Columns | New version | Historical unversioned exports |
| --- | --- | --- |
| `family_id`, `hypothesis_id` | Full family and exact retained hypothesis; align with generic q-values and size | Six-trait family and focal trait slug; do not align with generic full-family q-values and size |
| `full_family_id`, `full_hypothesis_id` | Explicit aliases of the generic full identity | Only `full_family_id` was recorded; confirm exact full membership using `feature_id` and the saved declaration |
| `q_by`, `q_bh`, `family_size` | Full-family correction values | Full-family correction values |
| `family_id_six_trait`, `hypothesis_id_six_trait` | Exact six-trait correction identity | Not recorded separately; the generic IDs have this meaning |
| `q_by_six_trait`, `q_bh_six_trait`, `family_size_six_trait` | Six-trait correction values | Six-trait correction values |

`full_family_membership=member` identifies a retained full-family member even
when it is untestable and its reported p/q values are unavailable. A focal
feature absent from the retained vocabulary instead has
`full_family_membership=absent_from_retained_vocabulary`: `hypothesis_id` and
`full_hypothesis_id` are null in the export mapping and blank in TSV. Its full
family ID and size remain context, not a claim that it is a member. The fallback
biological `feature_id` does not manufacture a full-family hypothesis. The
separate six-trait identity still joins to its declared unavailable slot;
unavailable probabilities remain null/blank rather than zero.

For each full member, join `(family_id, hypothesis_id)` to the saved
family's member inventory and to its correction record's embedded
`family.family_id` and `results[].hypothesis_id`. Use
`(family_id_six_trait, hypothesis_id_six_trait)` for the separate six-trait
correction. Absent full members have explicit nonmembership instead of a
fabricated full correction join.

Published unversioned tables and their bound manifests are not rewritten. A
consumer migrating across versions must distinguish absent legacy markers from
the explicit new marker, use the corresponding identity mapping, and refuse
unknown versions or contradictory aliases. This producer change adds no generic
legacy/new table reader. New focal headers and identity values, version metadata,
README text and output hashes change; the full screen and composition output
contracts and scientific correction values remain unchanged. Any historical
table migration or regeneration requires a separate reviewed plan.

For each profile, deduplicate reported taxon IDs. Reconstruct raw counts by
summing their consumed matrix vectors, and require exact agreement with the core
profile counts. Binary trait-bearing-taxon counts and matched-taxon denominators
are read from the verified core `ora_profiles` product and independently
reconciled with the deduplicated reported taxa and consumed matrix. Three
representations are kept separate:

- `raw`: a trait's annotation count divided by total retained annotation counts.
- `presence`: a trait's binary taxon count divided by summed binary taxon-trait
  presences. A taxon bearing many traits contributes to many denominator entries;
  this removes repeated annotation depth but **retains annotation-breadth
  weighting**.
- `prevalence`: the same binary count `k` divided by the ORA product's
  `taxa_count`, the number `N` of distinct matrix-matched reported taxa, including
  matched taxa with no retained annotation. Taxa with no matrix match are audited
  separately and are not silently treated as trait-negative observations.

All three representations use the same complete
human 16S/WMS experiment pairs with positive retained totals on both sides.
The prevalence comparison deliberately retains this shared pair restriction;
it is not a rerun of every ORA test or a prevalence estimate across patients.
Unmatched taxa, matched taxa without retained annotations, unpaired profiles and
zero-total pairs are counted explicitly. Matrix-matched all-zero rows remain
distinct from lookup failures.

Primary profile depth is the owner's specified mean of
`log10(1 + maximum retained trait count)` over distinct matrix-matched reported
taxa, including matched all-zero rows. Subtract decreased from increased after
computing each profile mean. The labelled secondary covariate is
`log1p(total retained annotation counts / distinct matrix-matched reported taxa)`.
Both are recomputed under this exact vocabulary and matching rule; no historical
numerical result is imported.

Within each study, average the eligible paired fraction differences and the
paired depth differences. Every study then contributes one observation. This is
the mean of paired transformed depths, **not** the transform of a pooled or
study-mean depth. Fit separate models for the primary and secondary covariates:

`study_mean_fraction_difference = intercept + slope * study_mean_depth_difference`.

The intercept is the model's conditional difference at zero depth difference.
It is only reported as an inferential estimate when zero is within the observed
depth range, the design has sufficient rank, residual degrees of freedom exist,
and leverage permits covariance estimation. A constant depth difference of zero
reduces to an intercept-only model; a constant nonzero difference cannot identify
the reference intercept. Empty, insufficient, unsupported, singular and
zero-residual-variance cases remain explicit, never silent zeros.

Uncertainty uses the [MacKinnon–White HC3 covariance](https://www.statsmodels.org/stable/generated/statsmodels.stats.sandwich_covariance.cov_hc3.html)
at the study level: the sandwich covariance uses squared residuals divided by
`(1 - leverage)^2`. Two-sided t reference probabilities and 95% intervals use
`n_studies - model_rank` degrees of freedom; this is an approximation, not a
calibrated guarantee. The unadjusted mean uses the same intercept-only HC3
convention. The sign screen calls `trait_screen.run_screen` on the reconstructed
paired profile inputs, preserving its full-candidate BY/BH, effective-m BY
sensitivity, exact Clopper-Pearson interval and test-status fields. Its test is
the positive share among non-tied study mean differences against one half.
The shared `SIGN_TOLERANCE = 1e-12` classifies absolute study differences at or
below that tolerance as ties, as in the native directional-consistency screen.
Pairing alone does not remove between-list annotation-depth differences, and a
linear depth covariate does not establish that all curation bias has been removed.

The composition sensitivity calls `trait_composition.compute_log_ratio_shifts`
on the raw and binary-presence focal counts. It preserves the native positive
pair exclusion, equal-study aggregation, log2 ratio sign test, exact interval,
and eligible-study/pair counts. Prevalence uses the same within-profile `N` for
both traits, which cancels in a trait ratio, so a duplicate prevalence ratio is
not reported. The declared family retains all 15 unordered focal pairs with BY
and BH corrections. If a core lacks a focal feature, its affected pairs are
explicitly unavailable, never filled with zero. The native helper runs on the
available features; its corresponding family size and BY values are separately
labelled beside the authoritative complete 15-slot correction.

Run `just presence-depth-sensitivity --core-run RUN_JSON --output-dir OUTPUT`
after the core has verified. The recipe writes the full screen, six-trait comparison,
composition summary and study ledger,
profile coverage/depth audit, experiment eligibility ledger, compressed study
effects, declared/executed test families, analysis bases and provenance. All
artifacts are bound through separate count and ORA consumer records. Provenance
includes Git commit and dirty state, consumed product hashes, all canonical
`mtp` Python modules loaded in the process at provenance construction, the policy
artifact, and dependency lock. This inventory includes indirect verification and
aggregation imports; in a long-lived process it can also include modules loaded
by earlier operations. It describes source coverage, not line-level execution.
At entry to `run_sensitivity`, before core consumption or scientific computation,
the producer snapshots all canonical Python source files and the shared policy.
Before creating the output
directory, it compares the actually loaded source inventory with those original
bytes, including modules imported lazily during analysis. Changes or removal of
previously loaded source records abort publication. Unused source files are not
added to the reported loaded-module inventory. This guard rejects relevant source
or policy changes between run entry and publication (#953). It does not compare
against the bytes used at module import: edits after import but before run entry
are outside its detection window, especially in a long-lived Python process.
The `just` recipe starts a fresh process; its source inventory is not an attestation
of the interpreter's already-loaded bytecode (#964).
The lock captures package versions. Scope is read with the shared
`run_strain_rollup` API, and rollup vocabulary is
`ncbi_taxonomy_closure_rollup_<scope>` as in #828. The sensitivity supports `off`
and `phenotype`; a valid exploratory `full` core is explicitly rejected before
analysis because the multiplicity contract does not yet support that scope.
Outputs must
be outside the immutable core and in a new destination directory. The command
never reruns an arm or modifies the source bundle.
