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
BY as the primary correction and BH as a comparison. The focal table preserves
the full-family corrections and adds a separately named six-trait BH comparison
to match the Table 1 family. Missing focal features retain unavailable slots.
Family declarations and correction records use the authoritative LinkML
multiple-testing v1.1 models and `mtp.multiple_testing`; they are replay-verified
against their declared inputs before writing. The shared policy also reports BY
for the six-trait family, without replacing the full-screen primary call.
Candidates with all ties or zero residual variance are untestable, as are
unsupported or insufficient fits. Untestable candidates occupy their family
slots but receive no reported p-value or q-value; internal
correction uses p=1 for these slots. Comparisons across these sensitivity families
are exploratory and are not independent replication.

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
paired profile inputs, preserving its nominal full-family BY/BH, effective-m BY
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
includes Git commit and dirty state, consumed product hashes, the analysis and
family-adapter module hashes, native test/aggregation/correction/consumer helper
hashes, the shared tolerance source, policy artifact, and dependency lock. These
named source hashes describe the declared implementation coverage; the lock
captures package versions. Rollup vocabulary is
`ncbi_taxonomy_closure_rollup_<scope>` as in #828. Outputs must
be outside the immutable core and in a new destination directory. The command
never reruns an arm or modifies the source bundle.
