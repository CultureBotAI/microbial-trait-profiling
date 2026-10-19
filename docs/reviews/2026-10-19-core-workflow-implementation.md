# Execution-checked LinkML core workflow implementation

Implementation issue: [#629](https://github.com/CultureBotAI/microbial-trait-profiling/issues/629).
Implementation PR: [#634](https://github.com/CultureBotAI/microbial-trait-profiling/pull/634).
Branch: `629-execution-checked-linkml-workflow`, based on main. The original
working checkout and its existing analysis/figure changes were preserved;
implementation and verification use `/private/tmp/mtp-workflow-629`.

## Delivered behavior

- Tracked LinkML workflow, execution, artifact, validation and consumer records;
  generated models validate the records actually emitted.
- Tracked scientific models for source identity, flattened tables, catalog,
  count/ORA profile rows and coverage, reusing the existing Study and assertion
  model. Wide tables retain their existing physical representation.
- A single executable workflow specification drives seven stages and twenty
  products. Inputs, dependency paths, resolved settings and scientific units are
  explicit. Generated documentation is checked against the specification.
- Unique attempts, deterministic fingerprints, actual code/lock/input hashes,
  exact collection inventories, staged validation and atomic step publication.
  Failed attempts remain diagnostic records and cannot be consumed as success.
- Strict count contracts adapted from PR #224 without replacing the current
  binary ORA builder. Cross-product checks recompute profile values, matched
  denominators, coverage, cohort membership and table references.
- Source-cell lineage resolves direct versus descendant assertions and exposes
  reused source taxa across reported taxonomic ranks. Aggregate coverage
  explicitly leaves direct/descendant origin unresolved.
- Native screen, paired composition, Table 1 plots and ORA accept verified core
  manifests. Their output bindings include coverage for the actual retained
  feature set. Historical imports remain explicitly unverified.
- Offline CI execution, fault-injection tests and typed numerical comparison
  against existing outputs.

Commands and scientific definitions: [workflow guide](../workflow.md).
Product inventory: [generated product table](../workflow-products.md).

## Scientific checks

The synthetic genus G has three annotated descendant species, including A.
Reported observations G, G, A become distinct reported IDs G, A: count=4,
ORA=2, matched denominator=2. Source A contributes through both G and A;
the lineage query reports this overlap. The fixture also checks human versus
mouse hosts, MediaDive-only annotations, matrix matches without annotations,
taxa absent from the KG, a true missing matrix match, and incomplete pairs.

The supplied production KG requires explicit reconstruction of a continued
header and 190,505 continued agent-metadata fields. The parser preserves both
agent columns and reconstructs 6,134,514 complete edge records, excluding zero.
Its ordered stream of 1,150,386 taxonomy/trait triples is identical to the
stream used by the native numerical builder. The saved audit helper is
[`scripts/adhoc/audit_raw_kg_structure_629.py`](../../scripts/adhoc/audit_raw_kg_structure_629.py).

BugSigDB release identity is verified from actual v1.3.1 bytes. The KG inputs
are identified by their actual hashes and resolved file locations; the workflow
does not retroactively certify the KG download or infer its release from a path.

The production ingest was independently reconciled to the raw snapshot with
[`audit_production_ingest_629.py`](../../scripts/adhoc/audit_production_ingest_629.py).
The complete audit is saved locally at
`/private/tmp/mtp-workflow-629/.tmp/production-ingest-audit.json`.

| Quantity | Raw snapshot | Emitted study YAML |
| --- | ---: | ---: |
| Studies | 1,246 | 1,245 |
| Experiments | 4,457 | 4,435 |
| Signatures | 7,425 | 7,365 |
| Reported taxon occurrences | 59,486 source taxon chains | 59,486 leaf-taxon records |

All emitted signatures preserve the ordered leaf-taxonomy-ID lists in the raw
snapshot. There are no ID/name-chain count mismatches and no placeholder or
blank taxon names. The source field `NCBI Taxonomy IDs` separates distinct taxa
with semicolons and lineage levels with pipes; `MetaPhlAn taxon names` separates
distinct names with commas and lineage levels with pipes. Selecting the final
member of each lineage retains the reported taxon's most specific level; it
does not count every ancestor as another reported observation.

The 59,074 informational log entries were classified by reason and checked
against the source records.
Of these, 59,066 record this lineage-to-leaf conversion; the remaining eight
record rejected DOI strings. The high event count therefore does not indicate
missing taxa or placeholder substitution. The 60 skipped signatures all have
missing/`NA` source taxonomy IDs and missing/`NA` taxon names. No signature with
a nonempty source taxon chain was excluded. Those source omissions account for
the removal of 22 experiments and the single study `bsdb:37577672`. They describe
unavailable source annotations, not observed biological or analytical zeros.

Eight study DOI values are dropped by the inherited parser, which removes
`https://doi.org/` or `http://doi.org/` prefixes and then requires `10.`. Five
source values use a recoverable-looking `doi:` prefix: studies `25073603`,
`31609493`, `36209079`, `38383483`, and `39277616`; some also contain a trailing
period. Three others contain malformed source strings: study `29758946` has
`http://dx.doi. 438 org/10.3233/JAD-180176`, while `32842982` and `37980506` have
`0.1186/…` prefixes. These are metadata exclusions requiring separate curation,
not evidence that the corresponding studies or taxon observations were lost.

The ignored-file-inclusive comparison finds exactly the same 1,245 study YAML
filenames in the original `data/studies` collection, with every file
byte-identical to the isolated production ingest. Thus the DOI exclusions and
source-empty signature exclusions are inherited behavior, not numerical changes
introduced here. The original `data/derived/ingest_log.json`, however, contains
only two `bsdb:99` fixture events and cannot serve as production provenance;
the repository's existing test-isolation fixture documents that historical
overwrite. The new run's own complete log is bound to its ingest step.

The published artifacts from production attempt
`run-78d44d233b3f43348f50693774f375b4` were independently audited with
[`audit_core_profile_coverage_629.py`](../../scripts/adhoc/audit_core_profile_coverage_629.py).
Its input manifest is
`/private/tmp/mtp-workflow-629/.tmp/core-production/run-78d44d233b3f43348f50693774f375b4/run.json`
(SHA-256 `1efa9620f5fe8217b6dab57717706fff948da11dae7d415bb1d6aba2cdafdcc4`);
the saved results are
`/private/tmp/mtp-workflow-629/.tmp/production-profile-audit.json`.
These figures describe that completed attempt; final code and provenance
verification are recorded separately below.

The taxon matrix has 5,829 reported-taxon rows and 2,650 feature columns, of
which 2,333 are MediaDive medium features and 317 are retained by the nonmedia
filter. Of the taxon rows, 2,255 have no annotations in any feature, 163 have
only medium annotations, and 3,411 have at least one nonmedia annotation.
Thus removing medium features changes annotation coverage; carrying the
unfiltered coverage into a filtered analysis would overstate it. The consumer
coverage sidecars address this explicitly.

| Quantity | Count profiles: all hosts | ORA profiles: human hosts |
| --- | ---: | ---: |
| Profiles before pairing or zero-total exclusion | 7,329 | 5,829 |
| Profiles zero across all features | 227 | 195 |
| Additional profiles containing only medium annotations | 15 | 11 |
| Profiles zero after removing medium features | 242 | 206 |
| Reported distinct-taxon memberships, summed across profiles | 58,974 | 48,015 |
| Matrix-matched memberships without any feature annotation | 5,606 | 4,595 |
| Matrix-unmatched memberships | 0 | 0 |

Membership totals deduplicate taxon IDs within each profile but reuse taxa
across profiles; they are not numbers of independent organisms or unique
taxa. The 7,365 emitted signatures form 7,329 count profiles because multiple
signatures with the same study, experiment and direction are grouped. The
absence of matrix-unmatched IDs follows from the core builder emitting a row
for every reported taxon, including KG-absent taxa. It does not establish
complete KG annotation coverage. In particular, all-zero rows remain eligible
for the ORA matched-taxon denominator. No annotation is not evidence of
biological absence.

Core coverage contains 5,494 fully annotated, 1,608 partially annotated and
227 wholly unannotated count profiles; the corresponding human ORA counts are
4,400, 1,234 and 195. Here “fully annotated” means every matched reported taxon
has at least one feature annotation in the declared feature set, not that all
traits have been measured. All aggregate records correctly retain
`lineage_resolution=unresolved`; genus or family matrix rows cannot be
presented as verified direct genus or family assertions from this summary.

A real-data check in human study `bsdb:27115405`, experiment 1, illustrates
the different units. In the decreased set, `coccus shaped METPO:1000668` has
count-profile value 12, while its ORA profile value is 2 among 4 matched
reported taxa. `sphere shaped METPO:1000683` has count 2 and ORA count 1 in
that same set. The inherited annotation counts may exceed the matched-taxon
denominator; the ORA values do not. These are separate feature identities
and annotation aggregations, not measurements of cell morphology frequency
or microbial abundance in the sample.

## Verification results

The saved numerical comparison for production attempt
`run-78d44d233b3f43348f50693774f375b4` is
`/private/tmp/mtp-workflow-629/.tmp/validation/production-comparison.json`.
It compares actual table values with the original checkout's existing outputs;
the historical artifacts remain `imported_unverified` because matching values
do not establish their producing inputs or execution history.

| Product | Compared numeric cells | Changed cells | Dimensions and scope |
| --- | ---: | ---: | --- |
| Taxon matrix | 15,446,850 | 0 | Identical 5,829 rows and 2,650 features; also byte-identical |
| ORA profiles | 15,452,679 | 0 | Identical 5,829 human profiles, 2,650 features and matched-taxon denominator; also byte-identical |
| Count profiles, shared rows and columns | 15,196,203 | 0 | 5,829 shared profiles and 2,607 shared features; current product adds 1,500 rows and 43 features, removing none |

The count-profile expansion was independently joined to the published host
metadata and feature catalog with
[`audit_legacy_count_expansion_629.py`](../../scripts/adhoc/audit_legacy_count_expansion_629.py),
with complete results at
`/private/tmp/mtp-workflow-629/.tmp/validation/legacy-count-expansion-audit.json`.
Every added row is nonhuman, with no missing host labels: for example, 865
profiles are from mouse experiments, 165 from rat experiments and 111 from
dog experiments. The baseline's 5,829 keys are exactly the current product's
human keys. Thus the all-host count contract produces a broader cohort than
that historical file; users must preserve the intended host selection in
downstream analyses.

All 43 added columns are MediaDive medium features. Each is zero in every
human profile and nonzero in at least one nonhuman profile. The broader
feature vocabulary therefore introduces no nonzero cells into the human
cohort and does not add any nonmedia feature. This establishes the observed
scope difference; the reason the historical producer omitted those columns
cannot be proven from its unverified output alone. The 43 columns are retained
in the current count product because its contract carries the matrix's full
ordered feature vocabulary.

Repository verification passed:

- `just check`: 944 tests passed, 4 skipped, 1 deselected; generated models,
  typing, lint, LinkML graph/schema checks, Bandit, dependency declarations,
  secret scan and dependency vulnerability audit passed.
- `just precommit`: every hook passed. Git provenance uses GitPython's
  maintained native-Git interface; tests compare its results with Git for clean,
  dirty, detached and linked-worktree states and check that reading provenance
  does not rewrite the index.
- `just core-fixture`, `core-screen`, `core-ora`, `core-verify-consumer` and
  the selected-cell `core-lineage` recipe passed. The small core fixture has
  only one nonmedia trait; the production bundle supplies composition's full
  six-trait vocabulary.
- Regression tests execute native composition and Table 1 entry points with
  conflicting checkout KG metadata. Verified consumption uses the bound run's
  measured KG identity and leaves the unverified KG release unset.

The final complete checks are recorded in `.tmp/workflow-check-git-fix.log` and
`.tmp/precommit-git-fix.log`. The four skipped tests require production-derived
inputs absent from this isolated checkout; the
separate production execution and comparison above use actual local inputs.

A second full production execution completed as
`run-49a0c712ad3d4fc59d10804951fe511a`, at
`/private/tmp/mtp-workflow-629/.tmp/core-production-final/run-49a0c712ad3d4fc59d10804951fe511a/run.json`.
It records commit `3218b146a53db9d16ff3a85f33d36b8a55eeb3fc`; its matrix and both
profile files are byte-identical to the first execution. Its conservative
`git_dirty=true` value exposed the tracked-directory-symlink classification bug
in Dulwich, recorded in [#635](https://github.com/CultureBotAI/microbial-trait-profiling/issues/635).
The manifest remains unchanged. The final implementation uses native Git
semantics, and the regression cases also cover genuine untracked symlinks and
empty symlink targets. The diagnostic helper is retained under
[`scripts/adhoc/issue629_audit_git_state.py`](../../scripts/adhoc/issue629_audit_git_state.py).
This metadata correction changes neither the scientific operations nor their
validated numerical outputs.

The public production composition command also completed against the first
verified bundle. It read 7,329 profiles, excluded 1,500 nonhuman and 1,159
unpaired profiles, retained 4,670 profiles from 2,335 experiments, and used
317 nonmedia feature columns. Its statistics, provenance and checked core
consumption record are under
`/private/tmp/mtp-workflow-629/.tmp/consumer-smoke/production-composition/`.
The provenance records the run's actual KG hashes and leaves KG release
identity explicitly unverified.

CI and external review status are available on the implementation PR. No merge
was performed.

## Boundaries

This implementation standardizes explicit-I/O core execution, not operating
system file-access auditing or cryptographic producer authentication. Existing
legacy commands remain usable without claiming verified core provenance.
PR #575's deck builder can use the shared consumer-verification API; historical
slides require their producer chain to be rebuilt before acquiring verified
core lineage. No slides or scientific claims were silently relabeled.

Annotation representation is not measured physiological activity, cell count,
metabolite flux, or independent biological replication. Matrix matching alone
does not imply KG annotation evidence or a universal genus/family phenotype.
