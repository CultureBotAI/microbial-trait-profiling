# Presence and annotation-depth sensitivity (#772)

This additive analysis is declared before the coordinated rerun. It consumes a
verified core bundle, preserves its contrast population and already applied
direction orientation, and leaves the three arms' primary calls and plans
unchanged. It is an annotation sensitivity analysis, not a causal adjustment or
measurement of physiology.

The complete candidate family is the core catalog's ordered non-MediaDive
vocabulary, including candidates with zero observed signal. Its actual size and
fingerprint are recorded in an `AnalysisBasis`; no historical feature count is
hard-coded. Six focal traits are extracts of the full screen, not a newly
corrected six-test family. Each representation and test has a separately declared
full family, BY as the primary correction and BH as a comparison.
Candidates with all ties or zero residual variance are untestable, as are
unsupported or insufficient fits. Untestable candidates occupy their family
slots but receive no reported p-value or q-value; internal
correction uses p=1 for these slots. Comparisons across these sensitivity families
are exploratory and are not independent replication.

For each profile, deduplicate reported taxon IDs. Reconstruct raw counts by
summing their consumed matrix vectors, and require exact agreement with the core
profile counts. In the presence representation each matched reported taxon
contributes 0 or 1 per retained trait. Its profile fraction divides that trait's
presence count by the sum of presence counts over all retained traits, not by
patients or by the number of taxa. Both representations use the same complete
human 16S/WMS experiment pairs with positive retained totals on both sides.
Unmatched taxa, matched taxa without retained annotations, unpaired profiles and
zero-total pairs are counted explicitly. Matrix-matched all-zero rows remain
distinct from lookup failures.

Primary profile depth is total retained annotation counts divided by the number
of distinct matrix-matched reported taxa, including matched all-zero rows. Apply
`log1p` to that mean separately in each profile, then subtract decreased from
increased. The secondary covariate is the mean of `log10(1 + maximum retained
trait count)` over the same matched taxa, again differenced within experiment.
The latter resembles the historical exploratory metric but is recomputed under
this exact vocabulary and matching rule; no earlier numerical result is imported.

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
convention. The sign test separately tests the positive share among non-tied
study mean differences against one half, using an exact two-sided binomial test.
The shared `SIGN_TOLERANCE = 1e-12` classifies absolute study differences at or
below that tolerance as ties, as in the native directional-consistency screen.
Pairing alone does not remove between-list annotation-depth differences, and a
linear depth covariate does not establish that all curation bias has been removed.

Run `just presence-depth-sensitivity --core-run RUN_JSON --output-dir OUTPUT`
after the core has verified. The recipe writes the full screen, six-trait extract,
profile coverage/depth audit, experiment eligibility ledger, compressed study
effects, declared/executed test families, analysis bases and provenance. All
artifacts are bound to the verified core through a consumer record. Outputs must
be outside the immutable core and in a new destination directory. The command
never reruns an arm or modifies the source bundle.
