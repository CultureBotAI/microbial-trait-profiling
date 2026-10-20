# Binding imported assessments to an integrated analysis

Issue [#729](https://github.com/CultureBotAI/microbial-trait-profiling/issues/729)
prevents calibration or contrast-orientation results from a different input
corpus from supplying an assessment to an integrated evidence table.

The integration command checks all three SHA-256 fingerprints declared by each
supplied calibration/orientation assessment: the count profiles, experiment
table and raw BugSigDB dump. The check applies to both direct-input and
`--core-run` execution. Missing, malformed or mismatched fingerprints produce
an `unreadable` standing with the offending field, and no procedure verdict is
attached. A missing assessment remains `pending`; a compatible assessment keeps
its actual procedure-specific pass, fail or insufficient-budget findings.

Compatible historical imports retain their producer revision and version,
analysis source hash, calibration protocol identity or orientation review
metadata. The generated provenance includes the assessment's absolute path and
file hash, its claimed input fingerprints, the expected fingerprints and the
binding result under `assessment_status.<assessment>.source`. The assessment
sidecar also records the accepted source path/hash and binding status. These
identities describe the producer that made the imported file; they are not
replaced with the integration's current revision. Reading an assessment through
the low-level API without supplying expected input hashes is explicitly marked
`not_checked`; the integration CLI always supplies its actual input hashes.

This is an input-identity check, not independent authentication of a producer's
scientific claims, proof of protocol adequacy or clinical validation. A schema-
valid report can accurately record failed or incomplete calibration. In
particular, the current-core one-replicate directional report cannot certify
error rates. The full-budget protocol uses 1,400 replicates for each of its two
null scenarios and 400 for each of four injection scenarios, seed 182. Its
study-flip assumptions do not address dependence between studies or establish
independent replication of the selected six traits.

## Guild calibration eligibility

Issue [#730](https://github.com/CultureBotAI/microbial-trait-profiling/issues/730)
identified an additional procedure mismatch even with exactly matching inputs:
calibration-v2 retained zero-total experiment pairs for guild McNemar, whereas
the integrated production analysis excludes the entire pair when either retained
annotation total is zero. The current core contains 94 such pairs. Exclusion
retains 4,482 of 4,670 profiles and 939 of 959 studies; guild discordants change
from 420/290 to 420/293. These are comparisons of annotation-fraction directions,
not measurements of physiological activity.

Calibration-v3 applies the shared whole-pair exclusion inside every guild
replicate. It records `whole_pair_nonzero_retained_totals_v1` and the observed
support in its status artifact. All loaded studies remain in the resampling
universe so the deterministic flip assignments retain their historical order;
the actual tested support is recorded separately. Null direction swaps preserve
the unordered retained totals, so the same 94 pairs are excluded in every null
replicate. The seed, replicate budgets, injection design and acceptance criteria
are unchanged. The protocol amendment discloses the observed mismatch rather
than describing this repair as preregistered.

For a bound assessment, integration requires the recorded guild eligibility
policy before attaching its guild verdict. An older or mismatched policy yields
`not_applicable` for guild McNemar while preserving the historical findings in
the sidecar. The other procedure-specific results remain scoped separately.

The retained-output comparator, `just audit-assessment-refresh`, accepts explicit
`--before-calibration`, `--after-calibration`, `--before-integration`,
`--after-integration` directories and a fresh `--output` JSON path. It is scoped
to this v2-to-v3 repair: all non-guild calibration values and all integrated
columns except `calibration_status` must remain exact. Flip assignments, injection
summaries, feature mappings and family inventories must remain byte-identical.
Issue [#734](https://github.com/CultureBotAI/microbial-trait-profiling/issues/734)
also requires the intentionally changed guild outputs to retain unique result,
replicate and criterion identities, declared replicate counts and fixed policy.
Their eligibility-dependent tested counts and statistical outcomes may change.
This comparison complements the integration's source-binding checks, consumer
record verification and inspection of assessment/provenance sidecars; it does
not replace them or independently establish calibration validity.

## Influence applicability

The completed `190-native-v2` influence diagnostic is deliberately reported as
completed results without a procedure-level scientific verdict. It contains
lineage-specific exclusion results under fixed and recomputed backgrounds, and
those scopes must remain separate. Its historical native ORA uses propagated
annotation counts, positive-cell enrichment testing and a depletion fold cutoff
of 0.67. Current production ORA uses distinct reported-taxon presence, all
eligible profile/trait tests, and a depletion cutoff of 0.5. Historical native
composition inference also uses experiment-level tests, with the shared
study-weighted comparison labeled separately. Input equality alone does not make
these procedures interchangeable.

Consequently, completion of 25 influence scenarios cannot establish robustness
of every present procedure. The existing adapter preserves provenance and scoped
findings without manufacturing a pass. A current-procedure robustness assessment
requires an explicitly scoped evidential basis or a separately declared
production-parity influence protocol and rerun. Taxonomic removal remains an
influence diagnostic, not a causal intervention.

## Orientation follow-up

The two recovered/convalescent versus healthy-control classification disagreements
are already addressed by the unmerged implementation in
[PR #577](https://github.com/CultureBotAI/microbial-trait-profiling/pull/577),
tracked in [#574](https://github.com/CultureBotAI/microbial-trait-profiling/issues/574).
Issues [#597](https://github.com/CultureBotAI/microbial-trait-profiling/issues/597)
and [#630](https://github.com/CultureBotAI/microbial-trait-profiling/issues/630)
describe its required integration with the classifier changes already on main,
including a version bump and frozen-input comparison. Until that work lands and
orientation is regenerated, importing the matched historical assessment preserves
58/60 recorded agreement and two disagreements. The review was performed by a
coding agent and is explicitly unreviewed by a human; correcting these two labels
does not independently validate the remaining corpus.
