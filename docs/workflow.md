# Execution-checked core workflow

Issue [#629](https://github.com/CultureBotAI/microbial-trait-profiling/issues/629)
implements the raw BugSigDB snapshot to downstream core-product specification.
The authoritative graph is [`workflows/bugsigdb_core.yaml`](../workflows/bugsigdb_core.yaml).
The [product and stage tables](workflow-products.md) are generated from that file.
The Python operation registry provides implementations; it does not maintain a
second dependency graph. Just recipes expose the interface.

## Run and verify

Install the declared environment with `just sync`. Validate definitions with
`just workflow-check`, and run the offline synthetic fixture with
`just core-fixture`. Each attempt creates a new directory and prints its
`run.json` path. The fixture needs no downloads or production data.

Run existing local production inputs in a new output root:

```bash
just core-run workflows/bugsigdb_core.yaml \
  --input raw_dump=/absolute/path/data/raw/full_dump.csv \
  --input kg_edges=/absolute/path/data/kg-microbe/merged-kg_edges.tsv \
  --input kg_nodes=/absolute/path/data/kg-microbe/merged-kg_nodes.tsv \
  --parameter source_mode=pinned \
  --output-root /absolute/path/core-runs
just core-validate /absolute/path/core-runs/run-ID/run.json
just core-lineage /absolute/path/core-runs/run-ID/run.json ora_profiles
```

`source_mode=pinned` checks the actual BugSigDB bytes against the repository's
pinned release checksum and length. An arbitrary supplied file is not promoted
to a pinned release by a filename or user label. The source identity product
records actual consumed files and the dump's embedded timestamp. The external
raw snapshot is the workflow boundary: download execution predating this run
is not reconstructed. Existing download provenance remains a source declaration.
The separately downloaded live studies, experiments, and signatures CSVs are
not inputs to this core workflow. Only the full dump is ingested.

The supplied KG export used for the production check has a header continuation,
repeated agent-field names and
manual-agent metadata continued across physical lines. Its input contract
reconstructs these known records, assigns distinct names to repeated fields,
counts the repairs and rejects unfamiliar malformed rows. It preserves the
native builder's taxon/predicate/object tuples; it does not silently discard
short physical lines or rewrite the source files.

To query the fixture's overlapping genus and species annotations:

```bash
just core-lineage /absolute/path/run-ID/run.json count_profiles \
  --profile bsdb:1 bsdb:1-1 inc --feature METPO:1000608
```

The query verifies the run, reads the recorded KG inputs, and reconciles source
taxa with the stored matrix count. It distinguishes a direct assertion subject
from a descendant subject, returns supporting source edges, and reports source
taxa reused by multiple reported taxa. A higher-rank row with inherited positive
counts does not establish a universal genus or family phenotype.

## Contracts and biological units

[`schema/workflow.yaml`](../schema/workflow.yaml) defines workflow/product
definitions, parameter declarations, executions, artifacts, validations,
downstream bindings, and unverified historical imports. Its generated Pydantic
classes validate actual persisted records. Graph validation separately rejects
duplicate IDs, unresolved references, cycles, overlapping paths, duplicate
producers, incompatible bindings, and undeclared parameters.
[`schema/core_products.yaml`](../schema/core_products.yaml) defines scientific
row contracts alongside the existing Study and taxon/trait model. Wide TSVs
remain the efficient physical representation; validators stream their rows.
The legacy files named `*_long.tsv` are transposes, not tidy observation tables.

| Product | Counting and selection rule |
| --- | --- |
| Studies | Exactly the YAML members emitted by this ingest attempt, individually hashed and validated. Earlier YAML files cannot enter through a shared directory glob. |
| Flattened tables | Explicit study, experiment, signature and occurrence keys, parent references, and declared multiplicities. Taxon occurrences can repeat across signatures. |
| Feature catalog | Qualified feature identifiers retain predicates where needed. Culture-medium features are explicitly classified. Core products retain MediaDive. |
| Taxon-trait matrix | Each cell counts distinct source taxa at or below the reported row taxon with the annotation. Duplicate source edges do not create extra source taxa. |
| Count profiles | Sum matrix count vectors over distinct reported taxon IDs per study/experiment/direction. All hosts are retained by default; only 16S and shotgun (`WMS`) experiments are kept (#742). |
| ORA profiles | Count distinct reported IDs whose inherited matrix count is positive. `taxa_count` is the number of distinct reported IDs present in the matrix. Default host: `Homo sapiens`; only 16S and shotgun (`WMS`) experiments are kept (#742). |
| Coverage | Reported = matched + unmatched; matched = annotated + unannotated. Status and sorted unmatched/unannotated IDs accompany the retained-feature digest. |

Matrix matching and KG membership are different. The current matrix builder
emits every reported taxon, including a zero row when the KG provides no
annotation. Such a taxon is matrix-matched and remains in ORA's denominator.
No annotation is not evidence of biological absence. Aggregate coverage says
`lineage_resolution=unresolved`; the selected-cell query resolves direct versus
descendant evidence from the original KG instead of guessing from counts.

Count profiles can reuse the same source annotation through a reported ancestor
and its reported child. In the fixture, G has three annotated descendants,
including A, and the reported sequence is G, G, A. Deduplicating reported IDs
gives a count-profile value of 4 and an ORA value of 2 with N=2. These are
annotation representation units, not four organisms or independent observations.

The core's increased/decreased direction is the curated report for Group 1
relative to Group 0. It does not by itself establish disease-case versus
healthy-control orientation; downstream contrast classification must make that
determination explicitly.

Strict `count-input-v1` validation is adapted from issue #184 / PR #224 while
preserving the current binary ORA builder. Duplicate headers/keys, ragged rows,
blank counts, negative counts, and non-ASCII or noninteger count strings fail.
Legacy blank-as-zero handling is available only as an explicit reader mode;
the core workflow uses the complete-count contract.

Host filtering is an explicit parameter at the profile stage. Pairing,
zero-total exclusions, MediaDive removal, background construction and multiple
testing are downstream analysis choices; the raw core does not silently apply
them. A cohort override changes the recorded artifact cohort, and consumers
can reject a bundle that does not match their expected cohort.

## What execution records establish

Every attempt has a unique run ID. Its deterministic fingerprint binds the
workflow definition, actual external input content, resolved parameter values,
code inventory, lockfile and runtime identity. Time and output destination do
not define scientific identity. Source hashes include the working tree's actual
Python and schema bytes, so a Git revision alone cannot hide uncommitted edits.
Only declared environment variables can supply parameter values.

Collection identity is a canonical sorted inventory of relative member paths,
sizes and SHA-256 hashes. Adding, removing or changing a YAML member changes
identity. The run embeds its definition and product bindings; it does not infer
provenance by searching for similarly named files.

Each stage writes to a fresh run-scoped staging directory on the publication
filesystem. Declared outputs are checked before the stage bundle is renamed
into place. The runtime verifies inputs before and after operations, rejects
undeclared outputs, and checks for changes to code, specification or lockfile
during execution. A failed attempt persists diagnostics but cannot be opened
as a complete core bundle. Successful earlier stages remain available for
inspection; they are not mislabeled as a complete downstream run.

`core-validate` has three distinct responsibilities: schema/graph validation,
actual artifact and inventory verification, and scientific cross-product
validation. The last reconciles study membership and metadata with raw source
records and documented exclusions, flattened exports with the Study collection,
and matrix annotation counts and metadata with the consumed KG. It also
recomputes profile aggregates and denominators from the recorded reported-taxon
and matrix products and checks their catalogs, cohort membership, coverage,
keys and parent references. Source identity claims must agree with actual bytes.

These records establish consistency of explicit-I/O Python execution. They are
not an operating-system audit of every file open, a cryptographic signature by
a trusted producer, or proof that an upstream biological annotation is true.
Historical runs can be read after checkout code changes; their recorded code
inventory remains evidence of the bytes used then. External input files must
remain accessible for complete revalidation.

## Downstream consumption and migration

### Owner rerun commands

The standard owner reruns require a verified core manifest and a fresh output
directory. They never run the old profile-builder prerequisites or infer their
inputs from the checkout's `data/derived` files. Prepare one shared core bundle
and use the printed `run.json` path for all three arms:

```bash
just rerun-core /absolute/path/core-runs \
  /absolute/path/data/raw/full_dump.csv \
  /absolute/path/data/kg-microbe/merged-kg_edges.tsv \
  /absolute/path/data/kg-microbe/merged-kg_nodes.tsv

just rerun-ilyas /absolute/path/core-runs/run-ID/run.json /absolute/path/reruns/ilyas
just rerun-noel  /absolute/path/core-runs/run-ID/run.json /absolute/path/reruns/noel
just rerun-jed   /absolute/path/core-runs/run-ID/run.json /absolute/path/reruns/jed

just rerun-verify /absolute/path/reruns/ilyas/arm_run.json
just rerun-verify /absolute/path/reruns/noel/arm_run.json
just rerun-verify /absolute/path/reruns/jed/arm_run.json
```

`rerun-core` defaults to `source_mode=pinned` and verifies the actual BugSigDB
release bytes. The final optional argument can select another declared source
mode; an unverified file never becomes a verified release through this option.
Use a checkout containing the merged implementation and run `just sync` first.
The Noel native influence publication stage also enforces its existing clean-code
and pinned-KG protocol requirements; incompatible inputs fail with diagnostics.

| Command | Standard suite |
| --- | --- |
| `rerun-ilyas` | Enrichment, depletion, full rankings, background strata, ORA calibration, and standard Figures 1–5 with a freshly computed density audit. |
| `rerun-noel` | Directional screen, shared study-flip assignments, directional calibrations, contrast orientation, taxon/native influence, validation evidence, and integrated evidence. |
| `rerun-jed` | Paired composition, log-ratio sensitivity, calibration, density/coupling audits, ordinary heatmaps, guild maps, and Table 1 maps. |

Stage-specific replicate counts, seeds and statistical thresholds retain their
existing defaults. `--replicates N` and `--seed N` explicitly override supported
stochastic stages; `--no-figures` omits optional renderers. Options are captured
in the execution record. Fixed scientific protocols remain explicit in the
individual command arguments. Figure rendering requires the existing rendering
dependencies, including Chrome for Plotly/Kaleido; install it separately with
`just get-chrome` if needed.

Noel's shared-flip artifact retains its fixed 100 draws. The screen generates
the reproducible study-flip sequence at the requested seed and replicate count;
the empirical power component uses at most 400 replicates. The existing
`matches_shared_flip_artefact` flag describes the seed-182/100-draw protocol,
so a custom seed or more than 100 screen replicates does not receive that label.

The original `scripts/rerun_*_arm_*.sh` entrypoints forward to these recipes and
now require the same two paths. A no-argument invocation fails with usage instead
of silently choosing old derived tables. Historical PD, six-trait exploration,
PPTX replacement, and branch-only directional scripts are outside this standard
suite. Their previous command lists remain in Git history; they must not be
represented as rerun by an `arm_run.json` from this suite.

Each attempt records its complete command plan before execution, actual command
exit status, logs, input/output checksums, code/configuration inventory, runtime,
and the upstream core manifest checksum in `arm_run.json`. The tracked LinkML
`ArmRun` class validates these records. A failed command leaves a failed record
and stops later stages; a failed attempt cannot pass `rerun-verify`. Outputs go
to the selected fresh directory, with verified core files copied into `inputs/`.
The source bundle and existing analysis outputs are preserved.

With the default population, count profiles retain all hosts in the core. Legacy
commands whose numerical algorithms assume human-only input receive an execution-recorded human subset
derived from `experiments.tsv`; no pairing, zero-count, or feature filtering is
added by that subset step. ORA uses the distinct human ORA count product. Native
consumers retain their own documented cohort and eligibility filters.

Final verification reconciles upstream core bytes, the ordered executed plan,
stage file snapshots, and the exact final output inventory. A directory snapshot
may gain files in a later declared stage; its already-recorded members must stay
unchanged. This records explicit command execution and file contracts; it is
not an operating-system audit of every file read or evidence that an annotation
is biologically true. Noel's integrated-evidence output retains its existing
explicit placeholders for unavailable native cross-arm results; running that
command does not imply that it imported Ilyas's or Jed's separate result tables.

The cross-arm comparison imports them from the three verified arm runs:

```bash
just rerun-compare /absolute/path/core-runs/run-ID/run.json \
  /absolute/path/reruns/noel/arm_run.json /absolute/path/reruns/jed/arm_run.json \
  /absolute/path/reruns/ilyas/arm_run.json /absolute/path/reruns/comparison
just rerun-verify-comparison /absolute/path/reruns/comparison
```

Each `integrated_evidence` input is the file its producing step recorded in that
arm's `arm_run.json`: Ilyas's `enrichment` and `depletion`, Noel's
`directional_screen`, `empirical_calibration`, `contrast_orientation` and
`native_taxon_influence`, and Jed's `composition_sensitivity`.
`--require-complete-sources` rejects a native table or status input that is missing
or unreadable, or a status that records other input bytes. No step produces the #184
coverage status yet, so it is left out and recorded as not supplied.

The `ComparisonRun` calibration verdicts and its eligibility ledger are built before
the stage writes anything, so a bad input leaves no partial directory behind (#844).

- **Noel's verdict.** It is read from `calibration_status.json`. That file covers
  `fraction_sign_test` and `screened_fraction_sign_test`, and each overall status must
  agree with its criteria.
- **The screened procedure.** It is the protocol's negative control, standing in for
  the selection of the six traits on this corpus, so the Table 1 verdict records its
  failure as anti-conservative. While it does not fail, the criteria are too weak to
  be evidence, so the Table 1 null control is never recorded as passed.
- **Jed's budget.** It counts the replicates in `composition_calibration_replicates.tsv`
  and requires them to agree with `composition_calibration_summary.json`.
- **The ledger.** It lists the studies that the shared basis keeps after the zero-total
  pair exclusion, computed from the core's count profiles. After integrated evidence
  runs, its recorded exclusion must match the ledger, with studies before equal to
  studies after plus studies lost (#814, #815, #844).

### Contrast populations and disease-oriented inputs

Select a population once when preparing the core. Use fresh destinations for the
core and every arm; the canonical checkout tables are not overwritten:

```bash
just rerun-core /tmp/disease-core \
  /absolute/path/full_dump.csv \
  /absolute/path/merged-kg_edges.tsv \
  /absolute/path/merged-kg_nodes.tsv \
  --contrast-population disease_vs_healthy_tier1_2_subject_controls

# Substitute the run.json printed by rerun-core:
just rerun-jed /tmp/disease-core/run-ID/run.json /tmp/disease-jed
just rerun-noel /tmp/disease-core/run-ID/run.json /tmp/disease-noel
just rerun-ilyas /tmp/disease-core/run-ID/run.json /tmp/disease-ilyas
just rerun-verify /tmp/disease-jed/arm_run.json
just rerun-verify /tmp/disease-noel/arm_run.json
just rerun-verify /tmp/disease-ilyas/arm_run.json
```

`just core-run --contrast-population ...` exposes the same selection through the
general workflow CLI. The optional positional source mode still works before the
flag, for example `fixture --contrast-population ...`; omitting it keeps the
verified pinned-source requirement.

| Population | Selection and direction |
| --- | --- |
| `all_paired` (default) | Preserve the existing unfiltered core and original Group 1 versus Group 0 directions. Human pairing and other eligibility remain downstream rules; the profiles still keep 16S and shotgun experiments only (#742). The name does not assert that every ingested record is human or paired. |
| `disease_vs_healthy_tier1_2_subject_controls` | Paired human disease-versus-healthy classifications with subject controls and Tier 1 healthy/normal or Tier 2 condition-negative/screened comparisons. |
| `case_vs_control_subject` | Paired human disease-versus-healthy classifications with subject controls, including Tier 3 hospital/comorbid controls. |
| `disease_vs_healthy_oriented` | All paired human disease-versus-healthy classifications, including paired-tissue and within-disease comorbidity comparisons. |

The three non-default selections require both ingested signature directions and
a 16S or shotgun (`WMS`) sequencing type; an experiment that fails only the latter
is recorded with exclusion reason `sequencing_type` (#742).
They do not impose trait-total eligibility, a strong-evidence-only filter or a
treated-disease flag exclusion. Later analysis filters can reduce support further.
Tier 2 condition-negative and Tier 3 subject controls are not claims of general
health; the broad population is not a pure disease-versus-normal subject cohort.
Counts are measured from execution, never fixed at a historical 959-study total.

The tracked [subset graph](../workflows/bugsigdb_core_subset.yaml) retains the
original ingested `source_studies` and complete source exports, then publishes an
exact selected Study collection plus `contrast_selection.json`. The selector uses
the classifier's `read_experiment_records` and `classify` APIs (there is no
`classify_corpus` API at this revision), and applies `oriented_profiles` to the
selected direction identities. Reversed experiments swap increased/decreased
once. Original group names, definitions, sample sizes, signature IDs and repeated
taxon occurrences remain intact.

Canonical experiment, signature, reported-taxon, count-profile, ORA-profile and
coverage products all derive from that same selected collection. Their direction
is disease relative to the selected comparison. The trait matrix and feature
catalog still derive from the **full source population**, preserving the ORA
background and vocabulary. Subset statistics are recomputed; they are not a
post-filter of the historical discoveries.

The core manifest's LinkML class is `WorkflowRun`. Both it and `ArmRun` record
`contrast_population` and `direction_basis`; verification checks these against
the executed parameter, selection stage and upstream core. The selection audit
records inclusion/exclusion decisions, reversed IDs, source and output inventories,
raw-input fingerprint and classifier version/source fingerprint. Its validator
reconstructs the permitted projection and rejects missing, added or altered
records, taxa, metadata, directions and audit values.

Noel's non-default workflow route uses an assessment adapter that requires a
verified core manifest and consumes its already oriented profiles directly.
It retains original raw classifications in its manifest and binds the exact
profile, experiment, raw-input and selection products from that core;
it does not reverse those profiles again. Hand-review counts are limited to
reviewed experiments inside the analysed population. This is classification
evidence, not physiological validation.

Run from the intended classifier revision. Selection does not repair classifier
rules, and verification uses the recorded classifier contract against the executing
implementation. After classifier changes, regenerate the subset core rather than
relabel an older manifest. Track 2 owns those rule changes and the final production
subset rerun.

The native trait screen, paired composition analysis, guild and Table 1 trait plots,
validation-evidence inventory, and ORA command accept `--core-run /absolute/path/run.json`. They obtain input paths
from the verified run and write a `core_consumption.json` binding (ORA places
it beside the report with the `.core_consumption.json` suffix). Conflicting
explicit input paths fail. Count and ORA profiles are not interchangeable.
The binding contains the run-manifest hash, product/artifact references, cohort,
feature policy and hashes of the downstream outputs.

Core consumers reject output paths that collide with their binding or coverage
sidecars, or overlap the immutable core run and source files. Directory-based
consumers require a dedicated analysis output directory; existing symlink or
hardlink aliases cannot redirect writes into protected artifacts.

Each binding also includes a checked annotation-coverage TSV for the consumer's
actual retained feature vocabulary. Thus a MediaDive-only taxon is recorded as
matched without retained annotations when the analysis excludes MediaDive.
This coverage describes the input cohort before analysis-specific pairing,
host/body-site restrictions or zero-total exclusions; those exclusions remain
in the analysis's own statistics and provenance.

```bash
just core-screen /absolute/path/run.json /absolute/path/screen
just core-composition /absolute/path/run.json /absolute/path/composition
just core-ora /absolute/path/run.json /absolute/path/ora
just core-ora /absolute/path/run.json /absolute/path/ora-depletion lower
just core-strata /absolute/path/run.json /absolute/path/strata
just core-verify-consumer /absolute/path/screen/core_consumption.json
```

`core-ora` writes `report.txt`, the selected `results.tsv`, the complete
`all_results.tsv`, and the report's consumer/coverage sidecars. Both tails use
global BH with q ≤ 0.001; enrichment requires fold ≥ 1.5 and depletion requires
fold ≤ 0.5. The complete export records the actual profile and taxa-matrix
SHA-256 values required by core-bound integrated evidence. Use a fresh output
directory when regenerating historical exports; adding fingerprints to an old
file does not establish that it was executed with those inputs.

`core-strata` binds human ORA profiles, experiments and the taxa matrix from the
verified core. It writes pooled, feces and feces-excluded partitions, complete
and selected ORA results, rankings, complete-trait summaries, density diagnostics
and sensitivity summaries. Its `core_consumption.json` inventories all these
products and the coverage sidecar; verify it with `just core-verify-consumer`.
The coverage describes the core input cohort before anatomical-site and
paired-profile eligibility restrictions. Extra analysis arguments, such as
`--site saliva`, follow the output directory. Explicit input overrides must
resolve to the bound core paths.

Partitioning rejects conflicting experiment body sites and missing experiment
keys before writing partitions. A present experiment with an unknown body site
still belongs to `feces_excluded`, which means "not identified as feces", not
"confirmed non-fecal". Identical metadata rows may be collapsed. The feces and
feces-excluded partitions are disjoint and together retain the pooled profiles.
Direct-path commands remain available without verified core provenance. This
standalone binding complements the owner rerun orchestration in #658 / PR #660.

The shared `verify_consumer_record` API is the integration boundary for deck
builders such as PR #575: verify the input analysis binding before accepting
its results or figures. A legacy figure/deck manifest cannot acquire verified
core lineage merely because its filenames match new products. Rebuild its
analysis/figure producer using a verified core run first. Existing direct-path
analysis commands remain available for historical reproducibility; they do not
emit a verified core binding.

Preserve old provenance claims without endorsing them:

```bash
just core-import-legacy /absolute/path/old-provenance.yaml /absolute/path/legacy.json
just core-compare /absolute/path/run.json /absolute/path/old-derived \
  /absolute/path/comparison.json
```

The import uses the tracked `LegacyImportRecord` schema. It hashes the original
document and retains its claims as `imported_unverified`; nonexistent historical
commands and paths do not become execution records. No ignored figure schema is
needed to read or validate the new workflow's provenance.

The comparison matches typed numerical cells by IDs and qualified features,
reports added/removed rows and columns, and distinguishes MediaDive differences.
Even identical values do not establish a historical baseline's producer or
inputs. Comparison never overwrites baseline or current products.

## Checks and maintenance

`just check` includes generated-model consistency, schema lint, graph checks,
and the offline pytest fixtures. CI also executes `just core-fixture` and saves
its run evidence. Tests cover mixed ranks, duplicate observations/source edges,
MediaDive-only and zero-annotation taxa, missing matrix matches, incomplete
pairs, nonhuman hosts, malformed counts, stale collections, changed inputs,
failed publication and consumer incompatibility.

Change schema first, regenerate with `just gen-pydantic`, update the executable
specification, and run `just workflow-describe` to refresh the product table.
Change the workflow version when its scientific contract changes. Do not edit
run manifests to make an old result pass a new contract.
