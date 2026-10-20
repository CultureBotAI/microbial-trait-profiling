# Adversarial review of PR #634

Scope: execution-checked LinkML workflow and core-product implementation in
[PR #634](https://github.com/CultureBotAI/microbial-trait-profiling/pull/634),
initial head `ea04a5d38dd6ed65a62a56f4c38aece74ee0e46d`. Review fixes are
tested after integrating main `41ca900`. The review covers execution integrity,
scientific source fidelity, biological interpretation, and downstream bindings.

## Confirmed findings

| Priority | Issue | Failure before correction | Required behavior |
| --- | --- | --- | --- |
| P1 | [#651](https://github.com/CultureBotAI/microbial-trait-profiling/issues/651) | Schema-valid loss or alteration between raw dump, Study YAML and flattened tables can produce a completed, mutually consistent bundle. Dropping an entire valid fixture study still completes. | Reconcile study membership, signatures, taxon observations and exports against declared source records; account for normalization and exclusions; reject blank or placeholder taxon names. |
| P1 | [#652](https://github.com/CultureBotAI/microbial-trait-profiling/issues/652) | Changing a fixture genus matrix count from 3 to 4 completes all seven stages and inflates the count profile from 4 to 5. Profile recomputation uses the same corrupted matrix. | Reconcile matrix counts, vocabulary and labels with consumed KG and reported-taxon inputs before publication and during revalidation. |
| P2 | [#653](https://github.com/CultureBotAI/microbial-trait-profiling/issues/653) | Authoritative schema prose implies repeated reported-ID multiplicity and automatically interprets direction as cases versus controls. | Describe distinct reported-ID aggregation and Group 1 versus Group 0 direction; clinical orientation requires downstream classification. |
| P1 | [#654](https://github.com/CultureBotAI/microbial-trait-profiling/issues/654) | ORA exits successfully when its results path equals the generated coverage path, then overwrites results with coverage. Its own consumer verifier subsequently rejects duplicate output entries. | Reject canonical collisions among outputs, record, coverage and protected core paths before any write. |
| P2 | [#655](https://github.com/CultureBotAI/microbial-trait-profiling/issues/655) | Schema-valid source identity can claim verified status for fixture bytes, invent a release or change the snapshot header while passing validation. | Reconcile pin flag, status, release and snapshot header with actual consumed bytes, including exact release fields for pinned data. |

The source-boundary regressions initially produced 13 failures and one pass.
The matrix mutation regression initially failed because no exception was raised.
Four source-identity contradiction regressions initially failed for the same reason.
The public ORA collision reproduction is retained under
`.tmp/adversarial-review-634/ora-collision/` in the implementation worktree.

## Biological interpretation

Reported taxa are abundance-shift annotations from literature; KG annotations
describe organism traits, not measured physiology in the sampled community.
Repeated reported IDs are deduplicated within a profile. However, the same KG
source taxon can contribute through two distinct reported ancestor/descendant
IDs. For reported G, G, A, with three annotated descendants of G including A,
the count profile is 4; the ORA profile is 2 with matched-ID denominator 2.
Those are different annotation units, not organism or cell counts.

Core direction records the curated increased/decreased report in Group 1
relative to Group 0. It does not establish disease cases, healthy controls or a
causal effect. The schema descriptions and generated models are corrected
without changing numerical aggregation.

## Review boundaries

No confirmed defect was found in the default runtime's path constraints,
collection inventories, input-retarget detection, code/lock mutation checks,
failure rejection or atomic step publication. Its focused runtime and Git-state
suite passed 55 tests before the fixes.

Consumer records establish verified upstream bindings and output-byte identity.
They do not claim a complete downstream execution recipe or cryptographic
authorship. This remains an explicit limit, not evidence that a scientific
annotation is true. Historical baseline products remain unverified even when
their numerical values match the new execution.

## Validation and disposition

All five findings are corrected in the review-fix commit of PR #634. A second
independent read-only pass found two additional cases within #654: an outside
hard link to a step publication record and a dangling report link to a future
statistics file. Both failures were reproduced and fixed before the final
full-suite run. The dedicated collision suite now has 19 passing cases.

Final local validation:

- `just check`: 1,226 passed, 4 skipped, 1 deselected; generated models,
  schemas, typing, lint, Bandit, dependency declarations, secret scanning and
  dependency vulnerability auditing passed. Log:
  `.tmp/adversarial-review-634-check.log`.
- `just precommit`: full hook-stack evidence is saved at
  `.tmp/adversarial-review-634-precommit.log`; final CI and review status are
  linked from the PR.
- `just core-validate` passed for the actual production bundle
  `run-49a0c712ad3d4fc59d10804951fe511a`. This rechecked every declared artifact,
  raw-derived Study fields and exclusions, all seven exported tables, all
  matrix cells against KG inputs, profile aggregations and source claims.
  Log: `.tmp/adversarial-review-634-production-validate.log`.

The source-export validator streams large profile rows. An earlier targeted
audit was interrupted after observing excessive memory use in its initial
eager comparison; it is not counted as successful evidence. Its helper is
retained at
[`scripts/adhoc/issue651_audit_source_boundaries.py`](../../scripts/adhoc/issue651_audit_source_boundaries.py).
The complete public validation above passed after the streaming correction.

Source informational log messages remain structurally checked rather than
independently reconstructed; exclusion events are independently reconciled.
No numerical aggregation rule was changed by these validation fixes.

The original production comparison remains documented in
[`2026-10-19-core-workflow-implementation.md`](2026-10-19-core-workflow-implementation.md).
The user authorized merging the checked PR and deleting its branch. Copilot
requests remain affected by [#267](https://github.com/CultureBotAI/microbial-trait-profiling/issues/267);
the independent adversarial review and all reproduced findings are documented
here, and no unavailable automated review is represented as approval.
