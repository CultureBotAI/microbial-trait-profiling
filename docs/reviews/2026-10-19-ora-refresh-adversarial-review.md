# Adversarial review of the ORA provenance refresh

Date: 2026-10-19. Reviewed PR:
[#722](https://github.com/CultureBotAI/microbial-trait-profiling/pull/722),
initial head `fa5d234`. Parent issue:
[#720](https://github.com/CultureBotAI/microbial-trait-profiling/issues/720).

The review separately examined the streaming audit and its tests, the standard
ORA recipes and consumer binding, and the actual production/biological evidence.
Current `main` was incorporated before final validation. The production inputs
and regenerated analysis outputs remain the ones recorded in the
[execution report](2026-10-19-ora-provenance-refresh.md).

## Finding and resolution

**P2 — matching malformed complete-result tables could pass the audit.**
Tracked in [#724](https://github.com/CultureBotAI/microbial-trait-profiling/issues/724).
The original audit checked table equality, basic columns, row widths,
significance flags and aggregate counts, but did not enforce the complete-row
scientific contract. Identical changes in both historical/current files could
therefore conceal invalid probabilities, counts, identifiers or missing
complete-export columns. Twenty paired-mutation regressions reproduced this:
the original helper returned PASS for each malformed pair.

The repair keeps streaming validation and uses the shared ORA column constants.
It checks complete-result fields, numeric domains, nonempty identifiers,
direction, background/profile-count consistency, eligibility, correction scope
and significance decisions against the declared p/q and fold thresholds. Error
tests assert the intended diagnostic, so an unrelated failure cannot make the
regression appear fixed. Actual-producer positive cases cover both tails,
global/profile correction, disabled FDR and zero-discovery/zero-eligible output.
Positive-infinite fold cutoffs retain their producer meaning; row statistics
and probabilities must remain finite.

## Production and biological evidence

The evidence review rehashed all five complete ORA exports and the final
integrated table, assessments and provenance against the committed audit and
consumer record. The recorded input hashes and the local report/audit copies
also agree. Logs record successful independent verification of enrichment,
depletion and integrated consumer records.

The regenerated results retain 66,240 pooled enrichment discoveries and 23
pooled depletion discoveries. The evidence table retains its 1,912 records,
including 40 method results. The 1,780 assessment records preserve 25 calibration
details, 19 orientation details, 1,701 influence details and 35 procedure-status
rows. The one-replicate calibration FAIL/insufficient-budget findings, 58/60
orientation agreement with two unresolved disagreements, and absence of a
procedure-level influence verdict remain explicit.

No new workflow or biological-interpretation blocker was found. The findings
remain contrasts in organism-level trait annotations among reported taxa;
they do not establish physiological activity, disease causation or independent
replication across arms. Previously documented limits remain: historical inputs
under `/private/tmp` require retention, stratified exports use direct-path
execution, and prior calibration/orientation assessments were imported rather
than rerun.

## Verification and disposition

The strengthened audit is rerun against the same five-pair manifest; its output
is compared byte-for-byte with the committed audit snapshot. Regression results,
full `just check`, precommit and final-head CI status are recorded in PR #722,
along with the successful production audit comparison and review disposition.
The analysis producer and the production result files are unchanged by #724.

The Copilot request limitation is tracked in
[#267](https://github.com/CultureBotAI/microbial-trait-profiling/issues/267).
Independent review and CI results are reported as such; they are not represented
as a Copilot approval.
