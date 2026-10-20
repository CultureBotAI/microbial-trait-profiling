# Follow-up adversarial review of the ORA refresh

Date: 2026-10-19. Reviewed evidence:
[the original review](2026-10-19-ora-refresh-adversarial-review.md),
[execution report](2026-10-19-ora-provenance-refresh.md), actual core products,
five historical/current ORA pairs, and imported assessments.

Unless stated otherwise, follow-up artifact paths below are relative to
`data/derived/ora_provenance_refresh/issue720` in this repository.

The original review established that the regenerated exports preserved historical
results. It did not establish complete test populations, durable retention of all
dependencies, or applicability of every imported assessment. The follow-up found
six actionable defects. Review of the new comparison helper found one additional
defect during implementation; it is tracked separately in #734.

| Priority | Finding | Issue and repair |
| --- | --- | --- |
| P1 | A calibration or orientation assessment from another corpus could attach a verdict to current results. | [#729](https://github.com/CultureBotAI/microbial-trait-profiling/issues/729): compare profile, experiment and raw-source fingerprints against actual integration inputs; retain producer identity and reject missing or mismatched identity. |
| P1 | Guild calibration and integrated production used different zero-total eligibility rules despite identical input files. | [#730](https://github.com/CultureBotAI/microbial-trait-profiling/issues/730): apply whole-pair exclusion within each calibration replicate, amend the protocol transparently, and withhold historical guild verdicts lacking the matching eligibility contract. |
| P2 | Identical historical/current exports could conceal duplicate tests and missing or invented population members. | [#725](https://github.com/CultureBotAI/microbial-trait-profiling/issues/725): derive eligible profiles and traits from actual inputs and require every expected profile/trait test exactly once. |
| P2 | Conflicting experiment metadata could multiply profile rows; missing metadata could silently enter the non-feces partition. | [#728](https://github.com/CultureBotAI/microbial-trait-profiling/issues/728): require a complete many-to-one experiment join before publishing partitions. Present-but-unknown body sites retain their declared policy. |
| P2 | Standalone stratification did not bind its partitions and products to verified core execution. | [#727](https://github.com/CultureBotAI/microbial-trait-profiling/issues/727): add explicit core execution, conflict/alias preflight, and a consumer record covering partitions, analysis products and coverage. |
| P2 | Reverification depended on source inputs and other-arm artifacts under temporary directories. | [#726](https://github.com/CultureBotAI/microbial-trait-profiling/issues/726): retain original records and dependencies in a LinkML-described archive, with explicit relocation mappings and verification that never falls back to historical paths. |
| P2 | The new before/after audit ignored guild rows before checking their population, allowing lost or duplicated replicates to disappear from comparison. | [#734](https://github.com/CultureBotAI/microbial-trait-profiling/issues/734): preserve unique guild ledger identities and declared replicate counts while permitting the expected outcome changes. The finding preceded any production comparison. |

## Scientific significance

The population audit passes all five production pairs without changing the saved
audit result. It checks eligibility, the matrix background and complete test
membership; it does not recalculate observed profile counts, effect estimates,
hypergeometric probabilities or BH corrections. Preserving historical numerical
results and independently validating a statistical method remain separate claims.

The calibration discrepancy affects 94 experiment pairs and changes retained
study support from 959 to 939. For the pooled guild McNemar calculation, the
discordant study counts change from 420/290 to 420/293. The repaired calibration
adapter reproduces production's exact p-value on these inputs. The protocol
amendment keeps the original seed, budgets, randomization assignments and decision
criteria, and records what was known when the amendment was made. A corrected
full-budget result must be assessed on its own; the earlier guild verdict cannot
be relabeled as validation of the repaired population.

A valid zero-discovery ORA result is distinct from a zero-eligible result. Matched
reported taxa without a positive trait annotation still contribute to the declared
ORA taxon denominator; unmatched taxon IDs are not implicitly counted. Regression
coverage must preserve this distinction rather than treating missing trait signal
as absence of a matched reported taxon.

These analyses describe organism-level trait annotations among taxa reported
increased or decreased in BugSigDB. Agreement among arms sharing these inputs
does not establish independent replication, physiological activity or causation.

## Corrected full-budget calibration

The corrected `calibration-v3` run completed from clean execution revision
`58b2ce841a88725603624007ca43d12c4a541c9d`, with seed 182, both 1,400-replicate
null scenarios and all four 400-replicate injection scenarios. The actual status
record binds the profile, experiment, raw-dump, protocol and analysis-source
hashes, and records the whole-pair eligibility policy and retained support.

| Procedure | Null replicates with any adjusted p-value at most 0.05 | Exact one-sided 95% upper bound | Family criterion, tolerance 0.065 |
| --- | --- | --- | --- |
| Fraction sign test, six-trait BH family | 128 / 2,800 | 0.0527474 | Pass |
| Log-ratio sign test, 15-pair BY family | 38 / 2,800 | 0.0177512 | Pass |
| Guild McNemar, single test | 118 / 2,800 | 0.0489325 | Pass |
| Screened fraction negative control | 250 / 2,800 | 0.0986499 | Fail, as declared for this negative control |

The first three procedures meet their declared overall criteria. The screened
negative control remains failed and is not promoted to a validated procedure.
The historical v2 guild result was 119/2,800; its different eligibility population
is preserved as historical evidence, not relabeled as the corrected result.

The power criterion checks reporting and Monte Carlo interval width; it imposes
no minimum power. For example, the fraction test's reported power at 1.1-fold and
1.3-fold injected signals is approximately 0.018 and 0.087. A passing reporting
criterion therefore does not imply useful sensitivity to small signals. The
per-test level checks remain sanity checks under the specified direction-flip
null; they do not assess dependence between studies or validate the native ORA
and composition procedures.

The complete v3 products are under the refresh's `issue730-calibration/results/`;
the v2 comparison run is preserved separately under `issue729-calibration/results/`.

## Core-bound stratification results

The fresh production rerun preserves all 24 historical scientific products:
23 files are byte-identical, and the sensitivity text report differs only in its
output-directory references. The comparison is saved as
`issue727-strata-comparison.json` beside the new `issue727-strata/` directory.

Profile keys are unique within each partition. Feces and feces-excluded keys are
disjoint and their union equals the pooled input: 3,643 + 2,186 = 5,829 rows.
These are input profile rows before the enrichment pairing filter, not numbers
of eligible ORA tests, independent studies or biological replicates. The existing
policy for present-but-unknown body sites is preserved; absent or conflicting
experiment metadata now fails instead of silently changing a partition.

After the pairing filter, the corresponding eligible profile counts are
4,670 = 2,950 + 1,720. With 317 tested traits, the saved ORA tables contain
1,480,390 pooled tests, 935,150 feces tests and 545,240 feces-excluded tests.
The shared background contains 5,829 taxa, a different quantity that happens to
equal the pooled input-profile count.

## Refreshed integrated evidence

The retained `issue729-assessment-refresh-audit.json` comparison passes. Relative
to the original integration with smoke-budget calibration, the refreshed main
table retains 1,912 rows and the method table retains 40. Exactly 22
`calibration_status` values change in each table; every other column is identical,
including effects, p-values, adjustments, support, orientation, influence,
sensitivity and coverage statuses. Feature and family tables are byte-identical.

In the full v2/v3 calibration comparison, all non-guild rows are unchanged.
Guild row identities, declared replicate counts and fixed decision policy are
preserved while its statistical outcomes reflect the repaired eligibility rule.
The randomization-assignment and injection artifacts are byte-identical. These
are exact comparison claims, not a substitute for source binding, consumer
verification or inspection of assessment and provenance sidecars.

The new integration provenance records matching actual input fingerprints and
the matching guild eligibility contract for calibration v3. It retains the
historical orientation producer, 58/60 review agreement and two disagreements;
the imported influence results still carry no procedure-level robustness verdict.
No native ORA or composition calibration verdict is introduced by this refresh.

## Retained production evidence

The production archive command completed successfully and published
`data/derived/provenance_archives/issue726/historical-ora-refresh/archive.json`,
with archive ID `archive:398c1e82-f07a-4f30-af16-fbeb3a9a726b`. Creation verified the
retained core execution and all three original consumer chains before publishing
the manifest. Historical pooled exports, both historical and refreshed strata,
native composition inputs, the earlier Noel arm outputs and imported assessments
are retained as inventoried supplemental evidence. Original records and results
were preserved. The creation log is under the refresh's
`logs/issue726-archive.log`.

The archive closes the actual historical input and result paths needed for these
checks. Descriptive references to the historical calibration protocol and
orientation review table remain recoverable from their recorded source commit,
`99160fee6338ddb5e8f30ed1da06739fefff0599`; this is not an archived software
environment or a rerun of those historical assessments. Later follow-up analyses
have separate execution records and outputs.

The retained bytes live locally under the repository's gitignored `data/` tree.
These PRs publish the archive tooling and specification, not a release asset or
an off-machine backup of the production archive.

Both follow-up producers completed successfully and published new records:

| Producer | Consumer ID | Path relative to the refresh directory |
| --- | --- | --- |
| Core-bound stratification | `consumer:d9c1c379-7bf8-451a-82d7-c05ba790a3a2` | `issue727-strata/core_consumption.json` |
| Integration with calibration v3 | `consumer:7ee68dec-11ae-4b8c-9a54-147121ac9a6c` | `issue729-integrated/core_consumption.json` |

These consumer records bind core products, output hashes and coverage. Integration's
analysis revision and native/calibration/orientation input bindings are retained
in `integrated_evidence_provenance.json`, itself a hashed output. Reading that
provenance and the assessment sidecar is part of the documented review scope.

Independent `just core-verify-consumer` runs passed for both published records
against the original core. Their logs are `logs/issue727-verify-consumer.log` and
`issue729-integration-logs/verify-consumer.log`. Archive creation separately
revalidates both chains using the retained copies.

These records postdate the historical archive. Their separate archive target is
`data/derived/provenance_archives/issue726/review-follow-up/archive.json`.
Use `just evidence-verify` with that full manifest path to check the later core
and consumer chains through retained locations instead of the original temporary
directory. The standard creation command publishes this manifest only after
retained inventory, scientific core contracts and both consumer chains pass.
Its execution log is `logs/issue726-follow-up-archive.log` under the refresh
directory. The historical archive remains unchanged.

Explicit follow-up retention covers the new outputs, both calibration runs and
their identity/eligibility evidence, both comparison audits, baseline integration
and strata, native composition and ORA inputs, historical Noel assessments, and
selected protocol, review-table, source, lockfile and comparison-helper snapshots.
This retention inventory does not turn supplemental assessments into new validated
consumer executions or claim preservation of a complete software environment.

## Implementation and review record

| Repair | PR | Execution or implementation revision |
| --- | --- | --- |
| Complete ORA test-population audit | [#731](https://github.com/CultureBotAI/microbial-trait-profiling/pull/731) | `2a6e678a73b88f12bfd9cafa7d4b786485eb47e2` |
| Retained, independently verifiable evidence | [#732](https://github.com/CultureBotAI/microbial-trait-profiling/pull/732) | `08daa792c0b73e444e4f638d436912b6355fd96a` |
| Core-bound stratification and metadata validation | [#735](https://github.com/CultureBotAI/microbial-trait-profiling/pull/735) | `58359153755f70c8c356d00b95d25f514a6711d5` |
| Corpus-bound assessments, corrected guild calibration and refresh comparison | [#733](https://github.com/CultureBotAI/microbial-trait-profiling/pull/733) | Calibration/integration execution: `58b2ce841a88725603624007ca43d12c4a541c9d`; later commits retain the comparison helper and this report. |

Separate coding-agent adversarial reviews checked the population logic, archive trust
boundaries and actual dependency closure, metadata join and partition contract,
assessment applicability and calibration eligibility, and comparison helpers.
Regression tests cover deliberate missing, duplicate, mismatched and tampered
inputs. The #734 helper repair passed all 19 focused tests before it was used on
production results. Green software checks alone are not a biological validation.

The combined source at `56e782498dde99f19eb9e58ba37ff04d1340799b` passed the
complete standard `just check`: 1,842 tests passed, four production-data tests
were skipped, and one test was deselected. Schema, workflow and security gates
passed, followed by the full pre-commit stack. The merge commit preserves the
original calibration/integration and stratification execution revisions as
ancestors. Final CI and review disposition are recorded on the linked PRs.

Copilot review requests returned success but registered neither a request nor a
review, including a request after marking the archive PR ready. This remains the
existing service/configuration issue
[#267](https://github.com/CultureBotAI/microbial-trait-profiling/issues/267).
It is not reported as a Copilot approval. The user's explicit merge-and-delete
instruction, independent reviews and CI results are recorded separately.

## Limits that remain explicitly scoped

The two recovered/convalescent versus healthy-control classification disagreements
have an existing proposed repair in
[PR #577](https://github.com/CultureBotAI/microbial-trait-profiling/pull/577),
with integration requirements tracked in
[#597](https://github.com/CultureBotAI/microbial-trait-profiling/issues/597) and
[#630](https://github.com/CultureBotAI/microbial-trait-profiling/issues/630).
The imported orientation evidence retains 58/60 agreement and two unresolved
disagreements. Its coding-agent review is not human validation.

The historical influence diagnostic uses procedures that differ from current
production ORA and composition. Its completed scenarios therefore remain scoped
findings without a procedure-level robustness verdict. Input-hash agreement does
not resolve this methodological mismatch. See the
[assessment applicability specification](https://github.com/CultureBotAI/microbial-trait-profiling/blob/main/docs/integrated_assessment_provenance.md).

Directional calibration does not calibrate the native ORA or composition methods
or the aggregated profile-support rule, nor does study-level randomization assess
dependence between overlapping studies.
The six Table 1 traits were selected on this corpus;
correction within that selected family does not provide independent confirmation
or establish error control for the selection step. Each reported family is
corrected separately, so a union of discoveries across families has no combined
false-discovery-rate guarantee.
The broader battery remains tracked in
[#263](https://github.com/CultureBotAI/microbial-trait-profiling/issues/263).
Archive verification establishes retained byte identity and the declared core and
consumer contracts; it does not authenticate historical authorship or reconstruct
the historical software environment for reexecution.
