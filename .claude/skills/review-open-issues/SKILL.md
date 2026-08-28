---
name: review-open-issues
description: Sweep and prioritize bsdbng's complete open GitHub issue queue using current ingest, provenance, schema, trait-profile, and analysis-branch evidence. Fetches every open issue (never the 30-item default), builds a pipeline-stage dependency graph from BugSigDB vintage through the two trait-analysis arms, checks each claim against code and artifacts on the branch it is about (not just main), groups root-cause families, and assigns a priority tier plus a separate rerun-cost class. Use for full backlog triage or deciding what is genuinely urgent; do not use as permission to close issues, regenerate data, or implement fixes.
metadata:
  category: workflow
  requires_database: false
  requires_internet: true
  version: 1.0.0
---

# Review and prioritize open issues

Adapted from the MicroGrowLink `review-open-issues` skill (v1.1.0).

Produce a complete, dependency-aware triage of bsdbng's open issues. The
issue queue, the analysis branches, and the committed artifacts on `main` are
different surfaces: sweep the queue itself, then test every claim against the
code and data it is actually about.

This is a read-only review by default. It does not implement fixes, regenerate
`data/`, close or edit issues, change labels, delete branches, or maintain a
tracker unless the user separately authorizes that exact mutation.

**When to use**: the user asks to review, triage, or prioritize issues or the
backlog; asks what is genuinely urgent; asks whether a collaborator has
addressed their issues; or a review pass has just filed a batch of issues that
need sorting.

**When NOT to use**: picking the next unit of work to implement, or acting on a
single known issue. This skill produces a ranking, not a fix. It is expensive
enough that it should not run on every "what's next" question.

## Sources of truth

Use these before relying on an issue title or an old planning note:

- `AGENTS.md` — required issue/branch/PR workflow, branch-deletion and
  merge rules, the no-workaround rule, and the data-fidelity rules (placeholder
  taxa, silent drops, column traceability);
- `README.md` — on `main` the only provenance statement is under "File
  types" (raw CSVs are generated, not committed). A "BugSigDB data provenance"
  section describing the Zenodo release pin (and which wiki CSVs remain
  unpinned) lands with #115 / the v1.3.1 pin; until that merges, read the
  branch, not `main`;
- `src/bsdbng/download.py` — the pinned release constants and URLs are
  defined here; this is the code-level source of truth for the vintage;
- `data/raw/download_provenance.json` — what a given checkout actually
  retrieved: resolved release, DOI, ETags, per-file `pinned` flag.
  **Gitignored, local only**, produced by `just download`; `git show` will
  never return it, so vintage claims must be measured on the machine that ran
  the download or read from `download.py`;
- `schema/bsdbng.yaml` — the LinkML schema; `src/bsdbng/datamodel/` is
  generated from it and is never the place to look for intent;
- `data/derived/ingest_log.json` — every drop and transformation assumption
  the ingest made; high-volume entries are bugs, not noise. **Gitignored,
  local only**, produced by `just ingest`;
- `Justfile` — the only supported entry points (`download`, `ingest`,
  `study-tsvs`, `taxa-trait-matrix`, `study-trait-profiles`, `check`); a
  recipe that exists only on a branch is evidence about that branch;
- `src/bsdbng/` and `tests/` for actual behavior on `main`;
- the analysis branches' own methodology documents, which are **not on
  `main`**: `trait-visualization/figure_methods.md` and
  `trait-visualization/trait-profile-investigation/*.md` on the composition
  arm; `data/derived/enrichment_methodology_summary.txt` and the `# ` metadata
  header block of `data/derived/trait_enrichment_results.tsv` on the ORA arm
  (the `*_ranked_*.tsv` exports carry **no** header; trace them back through
  `scripts/export_trait_enrichment_rankings.py`).

Treat issue bodies and titles as claims, not current status. Read comments:
this repository records corrections, withdrawals, and narrowed residual scope
there. A merged PR is evidence only after its code and acceptance criteria are
checked.

## What makes this repo different

**1. The repository has one queue under two names.** `turbomam/bsdbng` is a
GitHub transfer redirect to `CultureBotAI/bsdbng`. `AGENTS.md` still says
`--repo turbomam/bsdbng`; both resolve to the same issues and PRs. Confirm with
`gh repo view turbomam/bsdbng --json nameWithOwner` and report the canonical
name once, then stop worrying about it.

**2. Most analysis work lives on unmerged branches, not `main`.** The two
trait-analysis arms — the **composition / paired-contrast arm** ("Analysis A"
in issue titles; `scripts/create_study_trait_profile_heatmap.py`,
`scripts/create_trait_profile_infographic.py`, `trait-visualization/`) and the
**ORA/enrichment arm** ("Analysis B"; `scripts/trait_enrichment.py`,
`scripts/trait_depletion.py`, `scripts/export_trait_enrichment_rankings.py`) —
each sit on their author's branch on the `culturebot` remote. Verify the A/B
label against the script paths an issue cites; the letters are convention, not
metadata. An issue about either
arm is a claim about *that branch tip*. Reading `main`, or a stale local
checkout, answers a different question. Always:

```bash
git fetch culturebot --prune
git for-each-ref --sort=-committerdate \
  --format='%(committerdate:short) %(refname:short) %(authorname)' refs/remotes/culturebot
git show culturebot/<branch>:<path>          # read the file the issue is about
git log --format='%h %ad %s' --date=short <reviewed-tip>..culturebot/<branch>
```

Record the exact tip SHA you read against. "No new commits since the tip that
was reviewed" is a complete, citable answer to "has X addressed this?".

**3. Committed artifacts drift from the code that claims to produce them.**
Result TSVs, ranked tables, and figures are committed alongside scripts and
routinely predate the current script (older column names, excluded traits still
present, metadata headers that disagree with the methodology doc). Before
citing a number from a committed file, check its `# ` header lines and the last
commit that touched it (`git log -1 --format=%h -- <path>` on the branch) against
the script's current defaults. "Cannot be regenerated from the branch as-is" is
a real and recurring disposition.

**4. The two arms share one upstream and must stay aligned.** Both consume
`src/bsdbng/study_trait_profiles.py` -> `data/derived/study_trait_profiles_wide.tsv`.
An upstream filter (human-only, paired-experiment-only, MediaDive exclusion)
applied in one arm but not the other silently makes the arms non-comparable.
Issues of the form "arm A does X, arm B does not" belong at the shared upstream
stage, not in either arm.

## Workflow

### 1. Fetch the entire queue

Confirm the repository, current count, labels, and full queue. Never silently
accept `gh`'s default 30-item limit.

```bash
gh repo view --json nameWithOwner,url,defaultBranchRef
gh issue list --state open --limit 5000 --json number | jq length
gh issue list --state open --limit 5000 \
  --json number,title,body,comments,labels,createdAt,updatedAt,author
gh label list --limit 200
```

State the exact number reviewed and whether coverage was complete. Read every
issue body and its comments; for a long queue, inspect related groups in
parallel but preserve one disposition per issue.

### 2. Build the dependency graph before assigning rank

Place every issue at the earliest affected stage:

```text
BugSigDB release vintage (Zenodo full_dump.csv pinned; wiki CSVs live)
  -> download + download_provenance.json
  -> ingest: CSV -> per-study YAML, taxon normalization, ingest_log.json
  -> LinkML schema + generated Pydantic models, validation
  -> study TSVs (bugsigdb_taxa.tsv, experiments join)
  -> KG-Microbe download + taxa_trait_matrix.tsv
  -> study_trait_profiles_wide.tsv (shared by both arms)
  -> composition arm (A)     | ORA/enrichment arm (B)
  -> rankings, heatmaps, infographics, cross-arm evidence table
  -> release assets, README claims, manuscript/deck statements
```

An upstream identity, vintage, or filtering problem invalidates everything
below it. Recommend fixing or auditing that root problem before polishing a
figure or rerunning an arm. Group issues that share a root cause, but do not
hide the individual issue numbers.

For each issue, record when applicable:

- pipeline stage, and which branch (or `main`) the claim is about;
- affected BugSigDB vintage, KG-Microbe snapshot, and input/output file
  identities (path + branch + last-touching commit);
- schema/identifier assumptions (CURIE prefixes, taxon IDs vs names,
  placeholder names, EFO/condition labels);
- study / experiment / signature / profile / taxon counts and the filter
  that produced them (human-only? paired-only? MediaDive excluded?);
- which arm(s) consume the affected artifact and whether the arms still agree;
- prerequisites, blockers, duplicates, and superseding issues;
- cheapest decisive evidence and acceptance test;
- execution class: read-only audit; local regeneration
  (`just ingest`, `just study-trait-profiles`; `download` and `ingest` are
  wrapped in `just timed` and log to the gitignored `logs/pipeline.log`,
  `study-trait-profiles` is not timed — measure it, do not guess); KG-Microbe
  re-download
  (large, ask first); full arm rerun on a branch; release cut.

### 3. Check current reality and staleness

For each issue or group representative:

- Search exact issue references in history and open PRs:

  ```bash
  git log --all --oneline --perl-regexp --grep '#<N>\b'
  gh pr list --state all --search '<N>' --limit 100
  ```

  The word boundary is required: `#12` must not match `#120`. GitHub PR search
  is only a lead; open each candidate and verify that it actually resolves the
  issue. Also search commit messages on the analysis branches for the issue
  number and for phrases like "Addressed Issue" — authors reference issues in
  commit titles more often than in PR bodies here.

- Use `rg` (on a worktree or via `git show`) to confirm that named paths,
  functions, recipes, flags, and constants still exist and behave as
  described on the branch in question. Inspect tests as well as implementation.
- Compare acceptance criteria with the change. A commit message that says
  "implement per-profile FDR" but only adds an opt-in flag with the old
  default has not resolved an issue asking for per-profile FDR results. If
  only part is fixed, retain the issue with a narrowed residual.
- Distinguish an observation from its action issue. Prefer closing a fully
  recorded observation as superseded when a separate open issue owns the only
  remaining work.
- Verify artifacts by content and provenance, not filenames or prose. A
  ranked-trait table without its FDR scope, threshold, filter set, input
  vintage, and producing commit is not a canonical result.

### 4. Apply data-fidelity and analysis stop-the-line checks

Treat these as P0 when live or externally consequential:

- placeholder taxon names (`taxon_<ID>`), silently dropped rows, or a
  high-volume ingest-log entry that the pipeline is treating as noise;
- mixed input vintages: `full_dump.csv` from one release joined against
  wiki CSVs from another day, or `download_provenance.json` disagreeing with
  what `data/raw/` contains;
- a filter applied in one arm but not the other (host species, paired
  experiments, MediaDive traits) while cross-arm comparisons are being made;
- a statistic that is wrong in kind — wrong sampling unit, pseudo-replication,
  uncorrected extreme-value rollup, structural zeros reported as depletion —
  feeding a figure, README claim, deck, or release note;
- a committed result artifact that the branch's own code cannot regenerate,
  when that artifact is what an outward claim cites;
- a workaround commit (`# noqa`, `# type: ignore`, `# nosec`, hardcoded
  paths, `pip install`) that `AGENTS.md` says must become an issue instead;
- a provenance or tooling defect that would make an imminent release or
  large regeneration unusable while appearing valid.

Read `AGENTS.md` "Data fidelity rules" before applying these; that section is
the contract.

### 5. Assign priority and execution order

Use priority for consequence and a separate readiness/cost annotation for
ordering.

- **P0 — stop the line.** Live data corruption, mixed vintages, silently wrong
  scientific output, externally misleading claims, or a blocker that must be
  resolved before an already-planned release or expensive regeneration.
- **P1 — important and schedulable.** Correctness, reproducibility, provenance,
  or cross-arm consistency gaps in active analysis paths; missing guards for a
  likely workflow; test coverage that leaves a published number unverified.
- **P2 — low-risk or historical.** Documentation drift, refactors, theoretical
  edge cases, optional audits, tooling niceties, and work confined to
  superseded outputs without active spillover.
- **CLOSE/UPDATE.** Fixed, superseded, duplicate, no-longer-applicable, or
  title materially broader than the remaining work. Cite the exact
  commit/PR/code or comment that supports the disposition.

Calibrate P0 sparingly. Then order work within/across tiers using:

1. upstream unblockers (vintage, ingest, shared profile matrix) before
   arm-specific fixes;
2. shared-upstream alignment before either arm's statistics;
3. statistical validity before figure polish;
4. read-only and local-regeneration falsifiers before large downloads or
   full reruns;
5. golden-value tests on shared statistical helpers before more analysis
   built on them;
6. combine issues only when one patch or one measured run genuinely
   satisfies each issue's acceptance criteria.

Do not prioritize by age, by who filed it, or by a `P0` string in a stale
title alone.

### 6. Report

Return a compact report with:

1. coverage: canonical repository, timestamp, number reviewed, completeness,
   and the branch tips (SHA) read for each arm;
2. top 2–3 next actions, including why they unblock later work;
3. a dependency-ordered P0/P1/P2 table with issue number, current status,
   evidence, blockers, execution class, and next acceptance test;
4. CLOSE/UPDATE candidates with specific evidence;
5. per-collaborator summary when asked ("has X addressed their issues"):
   branch tip, new commits since last review, issues touched by those commits,
   issues untouched — stated flatly, no inference about intent;
6. unresolved evidence gaps and cross-repository ownership (KG-Microbe,
   BugSigDBExports);
7. a short sequence showing which costly work must wait.

Call out old issues explicitly rather than silently dropping them. Separate
measured findings, code inspection, inference, and proposed/untested work.

## Conventions this skill enforces

- **Full-queue coverage, not first-page sampling.** State exactly how many
  issues were reviewed and whether coverage was complete.
- **Evidence over vibes.** Every CLOSE/UPDATE/duplicate recommendation cites a
  specific commit, PR, artifact, or code location — never "this looks done."
- **Read the branch the issue is about.** A verdict measured on `main` for a
  branch-resident claim is not a verdict. Cite the SHA.
- **P0 is rare.** If more than ~10% of the queue lands P0, the calibration is
  wrong; recheck. A stale `P0:` string in a title is not evidence.
- **Titles are claims and they drift.** Issues get retitled mid-life while
  staying open. Re-read titles at report time rather than trusting the ones
  fetched at the start of the sweep.
- **The queue moves during the sweep.** Parallel pushes can resolve issues
  while triage is in progress. Re-check the open set and branch tips
  immediately before reporting, and say so if they changed.
- **Branches with commits ahead of `main` are never "stale".** Report what
  the commits contain; `AGENTS.md` forbids calling them stale or recommending
  deletion.

## Measurement discipline

The recurring failure is not misreading evidence, it is mismeasuring it.
Before citing any of the following, confirm how it was obtained:

- **The wrong ref.** `cat data/derived/x.tsv` reads your checkout, which may be
  `main`, a different branch, or a dirty worktree. For branch claims use
  `git show culturebot/<branch>:<path>` and cite the tip SHA.
- **`gh`'s silent 30-item default.** `gh issue list` and `gh pr list` truncate
  without warning. Pass `--limit` and print the count.
- **`gh issue list --json` omits `comments`** unless explicitly requested.
  Corrections and narrowed scope live in comments here; a body-only fetch
  systematically overstates what is open.
- **Exit codes through pipes.** `cmd | tail -3; echo $?` reports `tail`'s
  status. Use `cmd >out 2>err; echo $?`, or `${PIPESTATUS[0]}`.
- **Whitespace-splitting file lists.** Use `git status --porcelain -z` when
  paths may contain spaces.
- **Glob patterns tested by shape.** `git check-ignore --no-index <path>` tests
  what a `.gitignore` pattern does; a regex on the pattern tests what it looks
  like. Only the first is evidence.
- **Truncated tool output.** Several checkers elide long lines. Re-read the
  cited file at the cited line before acting on it.
- **Artifact headers vs. doc prose.** When `enrichment_methodology_summary.txt`
  and the `# Total-Tests-Performed:` header disagree, the header describes the
  file; the prose describes an intention. Cite the header.
- **Backticks in a double-quoted `-m`.** Write commit messages and reports
  containing shell examples via `-F <file>` or a quoted heredoc (`<<'EOF'`).

## Notes and limitations

- `gh pr list --search "<N>"` matches the number anywhere in indexed text and
  returns unrelated PRs. Treat every hit as a lead and open it before citing
  it. Likewise `git log --grep '#<N>'` needs the `\b` anchor, and `\b` is
  only a word boundary under `--perl-regexp`.
- An issue may be fully addressed in code while its acceptance criteria are
  not. Partial fixes keep the issue open with a narrowed residual; say which
  part is done and which is not.
- Evidence recovery is sometimes impossible: an intermediate result TSV may
  never have been committed. Say so and recommend superseding the issue with
  a regeneration request rather than leaving it open indefinitely.
- KG-Microbe trait semantics (METPO terms, trait propagation down the
  taxonomy, MediaDive media) are owned upstream. Flag, do not adjudicate.
- No @-mentions in issue comments or reports without explicit per-mention
  authorization (standing rule). Name collaborators in plain prose.

## Related

- `AGENTS.md` — required workflow (`gh issue develop`, PR with `Closes #N`,
  `just check`), branch and merge rules, data-fidelity rules.
- `src/bsdbng/download.py` and the README provenance section once #115
  merges — what is and is not pinned.
- `docs/bugsigdb-access-methods.md`, `docs/bugsigdb-live-counts.md` — dated
  measurements of upstream access paths and counts.
- Issues #131–#135 — the 2026-08-22 cross-arm statistical review of both
  analysis arms; the review text itself is not committed, the issues are the
  record.

## Mutation boundary

Do not close, comment on, relabel, retitle, or create issues during the
review. If the user later asks to act, present the exact issue numbers and
proposed mutation first. Apply closures one issue at a time with cited
evidence; never treat general approval as authorization for an unattended
bulk-close.

Do not regenerate `data/`, re-download KG-Microbe, or rerun an analysis arm as
part of triage. A recommended `just` recipe is a proposal, not permission to
run it.

Do not delete or force-update any branch, merge any PR, or open
cross-repository issues. Do not use `@` mentions without explicit
authorization.
