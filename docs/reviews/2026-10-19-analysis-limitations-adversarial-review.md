# Adversarial review of the analysis limitation fixes

Date: 2026-10-19. Scope: PRs [#553](https://github.com/CultureBotAI/microbial-trait-profiling/pull/553),
[#699](https://github.com/CultureBotAI/microbial-trait-profiling/pull/699), and
[#697](https://github.com/CultureBotAI/microbial-trait-profiling/pull/697).
The user requested adversarial review, issue creation, fixes, merge and branch deletion.
The original checkout's unrelated edits and existing analysis outputs are preserved.

## Findings

| Priority | Issue | Defect | Required behavior |
| --- | --- | --- | --- |
| P1 | [#705](https://github.com/CultureBotAI/microbial-trait-profiling/issues/705) | ORA exports from a different annotation background could pass a core join when profile identifiers were unchanged. | Complete exports must record actual profile and taxa-matrix hashes; core integration must require exact matches. |
| P2 | [#706](https://github.com/CultureBotAI/microbial-trait-profiling/issues/706) | Legitimate native study tests with missing primary p-values/intervals caused integration to fail. | Preserve explicit unavailable primary-test states, available effects and independent secondary results without inventing significance. |
| P2 | [#707](https://github.com/CultureBotAI/microbial-trait-profiling/issues/707) | Influence completion could be accepted with missing native outputs or compact comparison counts inconsistent with its execution manifest. | Require the producer's artifact inventory and reconcile available full/compact table counts before reporting results available. |
| P2 | [#708](https://github.com/CultureBotAI/microbial-trait-profiling/issues/708) | Standalone complete stratification export silently omitted unsupported traits when its profiles schema was absent or unreadable. | Require and strictly validate the input trait universe before writing a complete summary. |
| P2 | [#709](https://github.com/CultureBotAI/microbial-trait-profiling/issues/709) | A supplied depletion export was silently ignored when enrichment evidence was absent. | Require complete upper-tail evidence before importing auxiliary lower-tail summaries; fail before writing outputs. |

The ORA finding follows directly from the producer's absence of input hashes and
the importer's identifier-only comparison to the core. The native-test finding
is a producer/consumer mismatch: experiment contrasts can vary while their
within-study means are all zero. The influence acceptance gap was present in
existing passing fixtures with incomplete output inventories. The exporter gap
was reproduced against an actual analyzer fixture: five input traits, including
one unsupported trait, became four exported traits with zero reported untested
traits when the source profiles were omitted.

The integrated ORA representation uses upper-tail results as its primary source
and lower-tail discoveries as auxiliary summaries. A real lower-tail export on
its own now receives a clear error requiring enrichment evidence, instead of
being hashed into provenance while its results disappear into missing-source
placeholders. A regression exercises this path with the actual ORA analyzer.

## Statistical and biological review

The native study-level BH family, selected log-ratio BY family and full-trait
directional screen remain separate. Profile ORA hit summaries do not acquire
aggregate p-values. Stratification uses complete eligible denominators and
retains null traits; source background counts are not profile-count sums.
Influence execution remains distinct from a procedure-level scientific verdict.

Calibration's ten-replicate benchmark and its declared-spike truth limitation
remain as documented in the [implementation report](2026-10-19-analysis-limitations-followup.md).
This review does not turn that benchmark into evidence of mixed-signal FDR
control. The annotation contrasts do not measure in vivo physiology, microbial
cell abundance or metabolite flux, and do not establish causality.

## Verification and delivery

The fixes are delivered in the existing three PRs, with the additional issues
linked in their closing clauses. Regression coverage uses the actual producers
where source contracts are involved, including changed ORA inputs with unchanged
profile IDs, degenerate study contrasts, incomplete native/compact influence
evidence, and missing/malformed stratification input schemas.

The earlier bounded production integration predates #705. Its consumer record
verifies the recorded core and output bytes, but cannot retrospectively supply
the missing ORA input fingerprints. Those outputs remain historical evidence;
their hashes must not be relabeled as if produced by the corrected exporter.
Future core-bound integration requires newly exported ORA results carrying the
matching source identities. This review's fixture executions exercise that
corrected producer-to-consumer contract.

Final local checks, CI, exact commit inclusion and merge results are recorded
in the PRs. The merge sequence also checks the combined source tree. This report
records the findings and their resolution contracts; it does not claim a test or
remote action succeeded before its execution.

Copilot review registration remains unavailable under
[#267](https://github.com/CultureBotAI/microbial-trait-profiling/issues/267).
Independent adversarial review is recorded separately; absence of a Copilot
review is not an approval. The user's explicit merge/deletion request authorizes
those actions after the findings are addressed and checks pass.
