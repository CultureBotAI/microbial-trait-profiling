# Agent-assisted development

This project welcomes contributions made with coding agents (Codex, Copilot,
Claude Code, ChatGPT, or similar). The same rules apply to agent-generated
code as to hand-written code: it must pass `just check` before committing.

## Branches, deletions, and remote actions

- **Never delete a remote branch without explicit user approval.** Even if
  a branch looks stale, merged, or orphaned — ask first. A branch with
  commits ahead of main has unmerged work, period.
- **Never call a branch "stale" if it has commits ahead of main.** Report
  what those commits contain and whether the content reached main through
  another path (cherry-pick, separate PR). Let the user decide.
- **Never merge a PR without confirming all pushed commits are included.**
  If you pushed 3 commits but the PR only shows 1, the other 2 will be
  lost on merge. Check `gh api repos/.../pulls/N/commits` against what
  you pushed before merging or approving.
- **Questions are not requests.** When the user asks "why does this branch
  exist?" or "what's ahead?" — answer the question. Do not interpret it
  as permission to delete, merge, or take any action.
- **No destructive remote actions without per-action approval.** This
  includes `git push origin --delete`, `gh api -X DELETE .../refs/...`,
  closing issues, and merging PRs. Each action needs its own explicit OK.

## Do not accept workarounds

If your agent works around a missing tool, a broken path, or a failing check
instead of fixing the root cause:

1. Stop and note what the agent worked around
2. Open an issue describing the workaround and what should have worked
3. Fix the root cause, or ask for help

### Examples of workarounds that should become issues

- Agent runs `python -c '...'` instead of `just <recipe>` because `just`
  isn't installed
- Agent skips `just check` because a check fails
- Agent hardcodes a path because a config file is missing
- Agent adds `# noqa`, `# type: ignore`, or `# nosec` instead of fixing
  the code
- Agent installs a dependency that isn't in `pyproject.toml`
- Agent uses `pip install` instead of `uv sync`
- Agent creates a `.py` script in the repo root instead of putting it in
  `src/bsdbng/`

### Why this matters

A workaround that lets one person keep working creates a trap for the next
person. If the agent can't use the standard workflow, the standard workflow
is broken and we need to fix it — not paper over it.

## Required workflow for every change

Every change requires an issue, a linked branch, and a linked PR:

1. **Create or identify a GitHub issue** before writing any code.
2. **Create a branch linked to the issue** using `gh issue develop`:
   ```bash
   gh issue develop <number> --repo CultureBotAI/microbial-trait-profiling --checkout
   ```
   This formally links the branch to the issue in GitHub's Development
   sidebar. Do not create branches manually with `git checkout -b`.
3. **Push and open a PR** after the first commit. The PR body must
   include `Closes #<number>` to formally link the PR to the issue.
4. **Run `just check`** — all checks must pass before requesting review.
5. **Wait for CI and Copilot review.** Address all review comments.
6. **Do not merge your own PR** unless you are the repo owner.

No commits without an issue. No branches without a linked issue.
No PRs without a linked issue. No merging without CI green and
review addressed.

## What agents should do

- Run `just check` after making changes
- Use Justfile recipes (`just download`, `just ingest`, etc.) instead of
  running Python commands directly
- Keep branches short-lived. If a branch hasn't been updated in a week,
  either push progress or close the PR with a comment explaining why.

## Data fidelity rules

Every data pipeline step must be checked for silent data loss or corruption:

- **Look at actual output before declaring success.** Open a YAML file.
  Check that taxon names are real names, not placeholder codes. Check that
  fields you expect to be populated aren't all null. If 77% of your output
  has placeholder values, the pipeline is broken — not "mostly working."
- **High-volume log entries are bugs, not noise.** If the ingest log has
  87,000 entries saying "using placeholder name," that's not an edge case.
  That's the majority of the data failing silently. Investigate before
  moving on.
- **Verify encoding assumptions against real data.** Don't assume two
  columns use the same delimiter. Don't assume a column labeled "EFO ID"
  only contains EFO terms. Check with code: count unique values, check
  prefixes, compare counts across columns.
- **The ingest will fail on placeholder values.** If any taxon gets a
  placeholder name (`taxon_<ID>`), the pipeline raises an error. This is
  intentional. Fix the parser, don't suppress the error.
- **Every field in the output YAML must be traceable to a specific CSV
  column.** If you can't say which column a value came from, the mapping
  is undocumented and fragile.

## Biological interpretation and review

Every biological claim, figure caption, and summary text generated in this repository must adhere to the project-wide biological interpretation guide in [`methods/biological_interpretation.md`](methods/biological_interpretation.md) and trait semantics in [`docs/table1_trait_semantics_audit.md`](docs/table1_trait_semantics_audit.md):

- **Distinguish annotation from measurement**: BugSigDB signatures record reported taxon abundance shifts from literature; KG-Microbe records organism-level trait annotations from microbiological databases. We measure the alignment between reported taxa and trait annotations, not in vivo physiological activity, cell density, or metabolite flux.
- **Directional terminology**: Use precise terms corresponding to the specific pipeline quantity:
  - *Experiment/signature direction*: taxa reported increased or decreased in cases vs controls.
  - *Study fraction difference*: $\Delta = f_{\text{increased}} - f_{\text{decreased}}$ (paired within study).
  - *Log-ratio direction*: $\ln(\text{trait}_A / \text{trait}_B)$ paired shift between increased and decreased sets.
  - *ORA enrichment/depletion*: statistical over- or under-representation of trait-bearing taxa relative to background.
- **No speculative causal or clinical leaps**: Never assert that microbes "produce more SCFA", "cause inflammation", or "drive disease" based solely on trait-annotation enrichment or sign tests without direct physiological validation (see [`methods/validation_feasibility.md`](methods/validation_feasibility.md)).
- **Defensible wording**: State findings as "taxa annotated as [trait] were more frequently reported increased than decreased" or "enriched among increased taxa". Do not use uncalibrated terms like "dysbiosis signature" or "pathogenic shift".

## After a release

When a new release is tagged, check the release assets:

1. Check for new releases: `gh release list --repo CultureBotAI/microbial-trait-profiling --limit 1`
2. **Ask the user before downloading** — releases can be large
3. Download selected assets:
   ```bash
   mkdir -p release-assets
   gh release download <tag> --repo CultureBotAI/microbial-trait-profiling --dir release-assets \
     --pattern 'bsdbng-studies-*.tar.gz' \
     --pattern 'bsdbng-stats-*.txt' \
     --pattern 'bsdbng-ingest-log-*.json'
   ```
4. Unpack archives:
   ```bash
   mkdir -p studies
   tar xzf release-assets/bsdbng-studies-*.tar.gz -C studies/
   ```
5. Assess the output:
   - Read `release-assets/bsdbng-stats-*.txt` — compare counts against previous release
   - Spot-check a study YAML — are taxon names real or placeholders?
   - Check the ingest log — are skip counts reasonable?
   - Compare study count against BugSigDB live (https://bugsigdb.org/Main_Page)
6. Report any anomalies to the user before they share the release

## What to tell your agent

If your agent doesn't know about this project's conventions, tell it:

> This project uses `uv` for Python dependencies and `just` for commands.
> Run `just check` to verify changes. Do not add `noqa`, `type: ignore`,
> or `nosec` comments. Do not work around missing tools — tell me what's
> missing so I can install it.
