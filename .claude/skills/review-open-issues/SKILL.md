---
name: review-open-issues
description: Review and prioritize microbial-trait-profiling's complete open GitHub issue queue against current main, relevant PR heads, and data provenance. Use for full backlog triage or full-queue issue reviews; not for picking the next task or acting on a single known issue. Produces evidence-backed priorities and dependencies; reviewing alone does not authorize fixes, GitHub changes, downloads, or analysis reruns.
metadata:
  category: workflow
  requires_database: false
  requires_internet: true
  version: 2.1.0
---

# Review open issues

Adapted from MediaIngredientMech's
`.claude/skills/review-open-issues/SKILL.md`, retaining this repository's
branch-aware scientific review. Use the full queue, comments, measured evidence,
and upstream dependencies rather than a first-page sample or an old handoff.
Codex discovers this same directory through `.agents/skills/review-open-issues`
(a symlink; on a Windows checkout with `core.symlinks=false` it is a text file).

**When not to use:** picking the next unit of work to implement, or acting on a
single known issue. This skill produces a ranking, not a fix, and a full sweep
is expensive enough that it should not run on every "what's next" question.

## Establish the reviewed state

Run from the repository root. Read `AGENTS.md`. Record the review date, working
tree status, local HEAD, canonical repository, default-branch name and SHA, and
open PR heads/bases. Use `CultureBotAI/microbial-trait-profiling` explicitly in
GitHub commands. Former repository names in issue links may redirect here;
do not count them as separate queues.

Do not switch or reset a dirty checkout. Inspect committed files with
`git show <verified-sha>:<path>`, or use a separate worktree if execution is
needed and authorized. Discover remotes with `git remote -v`; do not assume a
remote named `culturebot` exists. Fetch relevant refs without pruning, or use
GitHub read APIs. If remote state is unavailable, report the review as
provisional rather than treating a cached ref as current main.

Keep three states separate: shipped on main, implemented on an unmerged PR or
branch, and present only in local files. A branch with commits ahead of main
has unmerged work; do not call it stale or infer permission to delete it.

## Fetch the complete queue and comments

Save API responses in a temporary directory so terminal truncation does not
become a hidden sampling limit. This recipe paginates the issue endpoint,
removes PR records, and obtains an independent open-issue count:

```bash
set -euo pipefail
review_dir=$(mktemp -d "${TMPDIR:-/tmp}/mtp-issue-review.XXXXXX")
gh api --paginate --slurp \
  'repos/CultureBotAI/microbial-trait-profiling/issues?state=open&per_page=100' \
  > "$review_dir/issue-pages.json"
jq '[.[][] | select(has("pull_request") | not)]' \
  "$review_dir/issue-pages.json" > "$review_dir/issues.json"
jq length "$review_dir/issues.json"
gh api graphql -f query='query {
  repository(owner:"CultureBotAI", name:"microbial-trait-profiling") {
    nameWithOwner
    issues(states:OPEN) { totalCount }
    defaultBranchRef { name target { oid } }
  }
}' > "$review_dir/repository.json"
gh api --paginate --slurp \
  'repos/CultureBotAI/microbial-trait-profiling/pulls?state=open&per_page=100' \
  > "$review_dir/pull-pages.json"
# Every issue comment in the repository, paged in bulk rather than one call per
# issue; group by issue_url below.
gh api --paginate --slurp \
  'repos/CultureBotAI/microbial-trait-profiling/issues/comments?per_page=100' \
  > "$review_dir/comment-pages.json"
jq '[.[][]] | group_by(.issue_url)
    | map({issue: (.[0].issue_url | split("/") | last | tonumber), comments: .})' \
  "$review_dir/comment-pages.json" > "$review_dir/comments-by-issue.json"
echo "$review_dir"
```

Each Bash call may start a fresh shell, so `$review_dir` does not persist
between calls: run the block as one script, then reuse the printed literal path.
`--slurp` needs gh 2.48 or later.

Inspect exit statuses and compare the saved issue count with `totalCount`.
Reconcile mismatches, including concurrent changes. A large `--limit` alone is
not proof of completeness, and the repository's `open_issues_count` includes PRs.

Read every issue body and its complete comment history from
`comments-by-issue.json`. Check that each issue's saved comment count equals
its `comments` field in `issues.json`; refetch any issue that disagrees.
(`gh issue list --json` omits comments unless they are explicitly requested.)

Corrections, withdrawn claims, and narrowed residual scope often live in
comments. Follow referenced parent issues, PR reviews, and closed issues when
they explain remaining work. Save one disposition per open issue even when
grouping related issues. A group representative does not establish that every
member has the same acceptance criteria or resolution.

## Ground each claim in this pipeline

```text
BugSigDB input vintage and download provenance
  -> ingest and taxon/condition normalization
  -> study YAML, schema validation, study/experiment/signature TSVs
  -> KG-Microbe snapshot and taxonomic trait propagation
  -> taxa-trait matrix and shared study/experiment/direction profiles
  -> composition | enrichment | study-level directional consistency
  -> paired log ratios, support joins, focus cuts and figures
  -> methods, biological claims and release artifacts
```

Use current files rather than freezing branch names, counts, or merge status
into the review:

| Evidence | Repository locations and interpretation |
| --- | --- |
| Authoritative index | `methods/README.md` lists the figures, recipes and methods notes on main; read it first and treat this table as orientation, not a complete inventory. `MANUSCRIPT_AUDIT.md` records audited manuscript claims. |
| Source contract and vintage | `README.md`, `src/bsdbng/download.py`, `src/bsdbng/ingest.py`, `data/raw/download_provenance.json`; distinguish the pinned dump from live auxiliary exports and inspect which files ingest actually consumes. |
| Ingest fidelity | `schema/bsdbng.yaml`, `src/bsdbng/ingest.py`, `data/studies/`, `data/derived/ingest_log.json` (gitignored, local only); generated Pydantic models implement the schema, not a separate source contract. |
| Trait counts and joins | `src/bsdbng/study_tsv.py`, `src/bsdbng/taxa_trait_matrix.py`, `src/bsdbng/study_trait_profiles.py`, `data/kg-microbe/`, `data/derived/`; verify IDs, directions, taxonomic roll-up, annotation coverage and join losses. |
| Statistical and display contracts | `methods/three_arm_statistics.md`, `methods/table1_paired_log_ratios.md`, `methods/table1_log_ratio_validation.md`, `methods/table1_trait_condition_body_site_maps.md`, `methods/guild_condition_body_site_heatmap.md`, `src/bsdbng/trait_guilds.py`, `src/bsdbng/trait_composition.py`, `src/bsdbng/stats.py`, `src/bsdbng/guild_maps.py`, `src/bsdbng/table1_trait_maps.py` (`just table1-compositional-stats`), and their tests. |
| Directional-consistency arm | `src/bsdbng/trait_screen.py` (`just trait-direction-screen`), `methods/trait_direction_screen.md`, and the pre-registered `methods/trait_screen_manifest.json`, whose dated amendments record changes made after results were known. |
| Branch-specific analyses | Resolve the current ref from the issue or PR. Examples change as branches move: enrichment work includes `scripts/trait_enrichment.py`, `scripts/trait_depletion.py`, ranking exports and `src/bsdbng/enrichment_stats.py`; calibration and provenance work includes `src/bsdbng/calibration.py` or a `calibration/` package, and `src/bsdbng/corpus_identity.py`. Verify existence at the cited SHA. |
| Execution and release evidence | `Justfile`, `tests/`, `.github/workflows/`, `logs/pipeline.log` (gitignored, local only), and available release assets; distinguish unit checks, integration checks, actual renders and scientific reruns. |

Historical issue labels use A for composition and B for enrichment. Confirm
the cited implementation rather than relying on those letters. Include the
study-level directional-consistency framework and log-ratio extension when
relevant. Do not assume either historical arm is still unmerged.

Before claiming absence, search ignored and hidden files with
`rg --no-ignore --hidden` or an equivalent exhaustive traversal of the relevant
location; state its scope and any exclusions. Use `git ls-files`,
`git status --short`, and `git check-ignore --no-index` to distinguish tracked,
untracked and ignored artifacts. The download-provenance sidecar is not
deliberately ignored by the current rules, but that does not mean a particular
copy is committed. Do not assign a producing SHA to an untracked file merely
because it sits beside that checkout.

## Verify resolution, provenance and dependencies

- Search exact issue references, for example
  `git log --all --oneline --perl-regexp --grep '#252\b'` (`\b` is a word
  boundary only under `--perl-regexp`). `gh pr list --search '<N>'` matches the
  number anywhere in indexed text and returns unrelated PRs. Search results,
  commit messages and merged PR labels are leads; inspect the actual diff and
  acceptance criteria. Check whether a claimed fix reached current main.
- Verify each issue on the ref it concerns, then separately state its main
  status. A fix in an open PR is pending integration, not a closure candidate
  based solely on its existence. Reverts, partial fixes and changes in defaults
  can leave residual work even after a PR merges.
- For artifact claims, record input hashes/vintages, producing code and command,
  filters, trait universe, denominator/background, sampling unit, tested family,
  thresholds, and eligible IDs where available. Last-touching commits and
  filenames are provenance leads, not proof that current code produced a file.
  Check actual metadata: a ranking export may omit the source TSV's header.
- Compare row/key sets, not just counts. Distinguish studies, experiments,
  signatures, profiles, taxa, and trait assertions. Before trusting a helper's
  output, inspect its arguments and confirm the loaded matrix is nonempty and
  has the expected traits; an empty matrix is not a decisive negative result.
- Inspect open PR dependencies and overlapping files. Record base/head SHAs,
  required merge order, unresolved review comments
  (`gh api repos/CultureBotAI/microbial-trait-profiling/pulls/<N>/reviews` and
  `.../pulls/<N>/comments`), CI state (`gh pr checks <N>`), and whether one
  change removes an invariant another depends on. Mergeability or green CI
  alone does not resolve scientific blockers or complete a missing review.

Use the cheapest decisive evidence first. Reading existing outputs is in scope;
regenerating data, downloading releases or KG inputs, installing renderers,
running external services, or rerunning analyses requires task authorization.
When execution is authorized, use the relevant `just` recipe and declared `uv`
environment. Inspect recipe dependencies first: `study-tsvs`, for example,
invokes ingest and is not a read-only inspection. Report missing tools rather
than bypassing the standard workflow. Preserve real command exit statuses:
`cmd | tail -3; echo $?` reports `tail`'s status, so use `set -o pipefail`,
`${PIPESTATUS[0]}`, or redirect to a file before checking `$?`.

## Scientific checks that change priority

Read `AGENTS.md` "Data fidelity rules" before applying these; that section is
the contract. Prioritize demonstrated damage at the earliest affected stage:

- Placeholder taxa, silent row drops, wrong ontology/identifier parsing,
  substantial unexplained ingest warnings, or input hashes that disagree with
  provenance. Inspect actual YAML and tables, not just successful exit codes.
- Incompatible populations or filters across compared analyses: human-only
  scope, complete increased/decreased pairs, MediaDive exclusions, retained
  traits, zero eligibility, and study weights. Live auxiliary CSVs alone do not
  prove mixed-vintage analysis if that path never consumes them.
- Wrong statistical unit, taxonomic assertion inflation, pseudo-replication,
  mismatched backgrounds, or an undeclared multiple-testing family. Shared-data
  agreement is not independent replication; native p/q-values from different
  hypotheses are not interchangeable evidence votes.
- Conflation of zero annotation, missing coverage, untested, nonsignificant,
  masked, or absent-from-display values. For log ratios verify numerator,
  reference, whole-pair positive-count eligibility, exclusions and correction
  family against the current implementation and methods.
- Focus panels that misreport why a condition disappears. Separate top-N
  selection, pooling into Other, minimum-study masking and blank-row removal.
  Check exact labels and set coverage: cuts may overlap and need not partition
  the corpus. A broad text match is not the requested biological category.
- Body-site significance presented as condition-specific significance,
  descriptive profile-support rules presented as calibrated cell tests, or
  annotation fractions presented as measured abundance, activity or causation.
- A workaround commit (`# noqa`, `# type: ignore`, `# nosec`, hardcoded paths,
  `pip install`) that `AGENTS.md` says must become an issue instead.

A proposed method in a document is not an implemented or validated result.
Distinguish a code-path risk from demonstrated damage in saved or published
outputs. For upstream KG-Microbe or BugSigDB defects, identify evidence, owner
and downstream impact without claiming the local repo fixes the upstream data.

## Rank and report

Keep priority, disposition and execution cost separate:

- **P0:** demonstrated silent corruption, invalid active scientific claims,
  security exposure, or a correctness blocker to an imminent deliverable.
  Name the affected output; do not promote every possible defect to P0. If
  more than about 10% of the queue lands P0, the calibration is wrong; recheck.
- **P1:** reproducible remaining defects or important reproducibility and
  scientific gaps that can be scheduled.
- **P2:** lower-impact documentation, maintenance, optional improvements or
  historical-only issues without current downstream impact.

Dispositions include active, partially resolved, fixed on main, implemented
but unmerged, duplicate/superseded, externally blocked, and unverified. Recommend
closure only with evidence addressing the issue's remaining acceptance criteria.
A moved/deleted function or an old issue date does not by itself prove closure.

For each issue record its pipeline stage, reviewed SHA(s), evidence, residual
work, dependencies, and a decisive acceptance check. Annotate the next action's
cost: read-only inspection, targeted tests, local regeneration, external/large
download, full analysis rerun, or release work. Order upstream unblockers before
dependent statistics and statistics before figure polish. Group common causes
without dropping issue numbers or silently treating distinct work as duplicates.

Before reporting, refresh the open set, titles, update timestamps and relevant
branch tips. Review new/changed issues and comments; identify issues closed
during the sweep. If the queue keeps moving or evidence is unavailable, name
the cutoff and coverage gap rather than claiming a complete current review.

Return the date, repository and SHAs, reviewed count versus full-queue count,
evidence limitations, a ranked table covering every issue, closure/update
recommendations with exact evidence, and the top two or three next actions.
Distinguish measured results, code inspection, inference and untested proposals.
When asked about a collaborator's progress, report commits and residual work
without inferring intent.

## Action boundary

Reviewing alone does not implement fixes, mutate GitHub, regenerate artifacts,
or launch expensive runs. Return the review without an approval detour.

- During the review, do not close, comment on, relabel, retitle, or create
  issues. Recommendations go in the report.
- No @-mentions in comments, reports or commit messages without explicit
  per-mention authorization. Name people in plain prose instead.
- Never treat general approval of the review, or of earlier actions, as
  authorization for an unattended bulk close.
- If the user has also authorized specific follow-up work, carry out exactly
  that scope. Closures, merges and remote branch deletion still require the
  repository's per-action approval, and the rules above still apply.

Do not silently turn triage into a tracker update or cross-repository issue
campaign.
