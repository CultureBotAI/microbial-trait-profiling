# Retaining execution evidence

Issue [#726](https://github.com/CultureBotAI/microbial-trait-profiling/issues/726)
adds a relocatable archive for already completed core executions and their
consumer records. This addresses the temporary-path dependency identified in
the [ORA refresh review](reviews/2026-10-19-ora-refresh-adversarial-review.md).

## Keep cited outputs out of temporary storage

Run outputs that a record cites must not live under a temporary directory, and must
be archived with `just evidence-archive`, then checked with `just evidence-verify`,
before the record merges. The system can empty its temporary directories at any
time, including on restart. Run F's verified run directories were lost this way
before anyone archived them ([#1040](https://github.com/CultureBotAI/microbial-trait-profiling/issues/1040)).

`just rerun-core`, `just core-run`, `just rerun-ilyas`, `just rerun-noel`,
`just rerun-jed`, `just rerun-compare`, `just rerun-multiplicity` and
`just evidence-archive` therefore refuse an output root under `/private/tmp`,
`/var/folders`, `/private/var/folders`, `$TMPDIR` or Python's
`tempfile.gettempdir()`
([#1049](https://github.com/CultureBotAI/microbial-trait-profiling/issues/1049)).
They check the path as written and with every symlink resolved. A durable-looking
link into temporary storage is refused, and so is a temporary path that links
somewhere durable, because the record keeps the path it was given. The refusal
happens before anything is written.

`/tmp` is not on the fixed list, because Bandit's B108 check rejects that literal
in the source
([#1051](https://github.com/CultureBotAI/microbial-trait-profiling/issues/1051)).
On macOS, `/tmp` links to `/private/tmp` and is refused. On Linux, `/tmp` is
refused only when it is the process temporary directory, which is the default.
If `$TMPDIR` points somewhere else, a Linux `/tmp` output root is **not** refused,
so do not use one.

Pass `--allow-temporary-output` only for a test or a throwaway run whose outputs no
record will cite. The offline `just core-fixture` run is synthetic and is not checked.

The check applies to an archive's destination, not to the runs it retains. A run
written to temporary storage can still be archived to a durable location while it
exists.

## Create and verify an archive

Use a fresh destination outside every retained input directory:

```bash
just evidence-archive /absolute/path/retained-evidence \
  --core-run /original/core/run.json \
  --consumer /original/enrichment/report.core_consumption.json \
  --consumer /original/integration/core_consumption.json \
  --retain /original/historical-results \
  --retain /original/other-arm-inputs
just evidence-verify /absolute/path/retained-evidence/archive.json
```

Every consumer's core run, core products, external raw/KG inputs, step records,
output files and coverage are retained automatically. A shared ancestor or
dependency is copied once. `--retain` explicitly adds historical comparisons,
other-arm inputs, assessment artifacts, logs or reports. Arbitrary paths inside
free-form JSON or report prose are not treated as dependency declarations;
include the corresponding evidence explicitly. Supplemental evidence is
inventoried but its scientific conclusions are not certified by copying it.

The original run, source identity and consumer records remain byte-identical.
The LinkML `EvidenceArchive` and `ArchivedLocation` records in
[`schema/workflow.yaml`](../schema/workflow.yaml) describe where each original
location is retained, with file and collection inventories. Archival time is
separate from historical execution time. The records do not assert a new
analysis execution, cryptographic authorship or independent biological evidence.

Creation checks for source changes during copying, inventories hidden files,
and verifies the retained core scientific contracts and consumer chains before
publishing `archive.json`. The destination must be fresh and disjoint from the
sources. An incomplete attempt has no published archive manifest. Failed
attempts remove only their own newly created destination.

Move the entire archive directory as needed, then use `just evidence-verify`
again. Resolution uses only the explicit archive mapping; it does not fall back
to the historical locations even if they still exist. The verifier checks the
core graph, inputs, scientific cross-product relationships, coverage and consumer
output fingerprints using retained bytes. Source-identity path claims are checked
against the original recorded core bindings, not silently replaced with archival
paths. Missing, altered, undeclared, symlinked and hard-linked members fail.

To locate a retained input for a new analysis or comparison:

```bash
just evidence-resolve /relocated/evidence/archive.json /original/input.tsv
```

This prints the retained path after checking the archive inventory. A new analysis
must record its actual inputs and execution normally. An archival path mapping
does not make old calibration, orientation or influence assessments applicable to
a different corpus or procedure; those checks remain the consumer's responsibility.

## Historical ORA refresh retention

The production retention destination for issue #720 is:

```text
/Users/marcin/Documents/VIMSS/ontology/KG-Hub/KG-Microbe/microbial-trait-profiling/data/derived/provenance_archives/issue726/historical-ora-refresh
```

Its manifest lists the three explicitly verified consumer chains: pooled
enrichment, pooled depletion and `integrated_with_assessments`. The core run,
original raw dump and external KG inputs are included through those chains.
Explicit supplemental retention includes the historical pooled ORA exports,
historical and refreshed strata, the native composition inputs, the earlier Noel
arm directory, and the refresh's input manifest, audit and two reports. Their
original paths and complete inventories are recorded in `archive.json`.

The creation log is `logs/issue726-archive.log` under
`data/derived/ora_provenance_refresh/issue720`. Creation publishes the manifest
only after all retained bytes and declared core/consumer chains verify. Run
`just evidence-verify` on the published manifest to independently repeat those
checks. Supplemental records may mention old script locations or other historical
executions: those claims remain descriptive unless explicitly included as verified
consumer/core dependencies. This archive preserves the original refresh; later
follow-up analyses retain separate execution records and outputs.
