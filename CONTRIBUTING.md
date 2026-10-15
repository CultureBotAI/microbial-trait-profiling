# Contributing to microbial-trait-profiling

## Setup

### macOS / Linux

```bash
# Install uv (Python environment manager)
curl -LsSf https://astral.sh/uv/install.sh | sh

# Install just (command runner)
# macOS:
brew install just
# Linux:
curl --proto '=https' --tlsv1.2 -sSf https://just.systems/install.sh | bash -s -- --to ~/.local/bin

# Clone and set up
git clone https://github.com/CultureBotAI/microbial-trait-profiling.git
cd microbial-trait-profiling
uv sync --all-extras --group dev --group qa
just install-hooks    # install pre-commit hooks
just check            # verify everything works
```

### Windows

If you have nothing installed, follow these steps in order.

#### Step 1: Install Git for Windows

Download and install from https://gitforwindows.org

This gives you `git` and **Git Bash** — a bash shell that `just` needs.
Use Git Bash (not PowerShell or cmd) for all commands below.

#### Step 2: Install uv

Open **Git Bash** and run:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Or if you prefer PowerShell:

```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

Close and reopen your terminal after installing so `uv` is on your PATH.

#### Step 3: Install just

Option 1 (winget, from PowerShell):
```powershell
winget install Casey.Just
```

Option 2 (pre-built binary): download from
https://just.systems/man/en/packages.html and put it on your PATH.

#### Step 4: Verify

Open **Git Bash** and run:

```bash
git --version
uv --version
just --version
```

All three should print version numbers. If any fails, it's not on your PATH.

#### Step 5: Clone and check

```bash
git clone https://github.com/CultureBotAI/microbial-trait-profiling.git
cd microbial-trait-profiling
uv sync --all-extras --group dev --group qa
just check
```

`just check` should pass. If it doesn't, open an issue.

#### Troubleshooting

- **`just` can't find bash**: make sure you're running in Git Bash, not
  PowerShell or cmd. The Justfile uses `bash` as its shell.
- **uv not found after install**: close and reopen your terminal.
- **Execution policy errors**: if PowerShell blocks install scripts, you may
  need to run `Set-ExecutionPolicy -Scope CurrentUser RemoteSigned` once.
  This only affects your user account.
- **Your coding agent works around a missing tool**: don't let it. See
  [AGENTS.md](AGENTS.md) for why and what to do instead.

## API keys

Some features (LLM embedding search via linkml-store) require an OpenAI API key.

```bash
cp .env.example .env
# Edit .env and add your key
```

The `.env` file is gitignored and will not be committed. Never commit API keys.

## Making changes

### 1. Pick an issue

All work starts from a [GitHub issue](https://github.com/CultureBotAI/microbial-trait-profiling/issues).
Don't start work without one. If you want to do something that doesn't have an
issue, open one first.

If BugSigDB is unreachable while reproducing pipeline steps, see the README section
"If `just download` is unreachable" for using the release raw-data tarball instead.

### 2. Create a branch

```bash
git checkout main
git pull
gh issue develop <number> --repo CultureBotAI/microbial-trait-profiling \
  --name issue-<number>-<short-slug> --checkout
# Example: gh issue develop 30 --repo CultureBotAI/microbial-trait-profiling --name issue-30-aggregate-tables --checkout
```

### 3. Write code

- Source code goes in `src/bsdbng/`
- Tests go in `tests/`
- Justfile recipes for new pipeline steps should use the `timed` wrapper

### 4. Run checks before committing

```bash
just check
```

This runs everything: ruff (lint + format), mypy (types), pytest (tests),
linkml lint/validate (schema), bandit (security), deptry (dependencies),
and pip-audit (vulnerabilities).

**All checks must pass.** If something fails, fix it. If you think the
check is wrong, open an issue — don't skip the check.

If you want to auto-fix formatting:

```bash
just format
```

### 5. Follow biological interpretation standards

If your contribution adds or modifies prose describing biological results, figure captions, or methods:
- Adhere to [`methods/biological_interpretation.md`](methods/biological_interpretation.md) for approved phrasing, evidence chains, and caption contracts.
- Consult [`docs/table1_trait_semantics_audit.md`](docs/table1_trait_semantics_audit.md) for trait definitions and supported interpretations.
- Maintain distinction between reported taxon annotation patterns and in vivo physiological measurements.

### 6. Commit and push

```bash
git add <your files>
git commit -m "Short description of what and why"
git push -u origin issue-<number>-<short-slug>
```

### 7. Open a PR

- Reference the issue: "Closes #30" in the PR body
- Wait for CI to pass
- Respond to Copilot review comments if any

## What not to do

- Don't commit generated files (`data/raw/*.csv`, `data/studies/*.yaml`,
  `data/derived/*`, `logs/`) — they're gitignored for a reason. The
  exceptions are the small summary tables a methods note cites and the figure
  outputs under `trait-visualization/`; see
  [Generated outputs on `main`](#generated-outputs-on-main).
- Don't edit `src/mtp/datamodel/mtp_pydantic.py` or
  `src/bsdbng/datamodel/bsdbng_pydantic.py` by hand — they're generated from
  the schema. Run `just gen-pydantic` after schema changes.
- Don't add `# noqa`, `# type: ignore`, or `# nosec` comments. If the
  check is wrong, fix the code or open an issue about the check.
- Don't add dependencies without discussing them first. This project uses
  a minimal dependency set intentionally.
- Don't use sentinel values (e.g., `NCBITaxon:0`, `"unknown"`, `-1`).
  If data is missing, use `null` or skip the record and log why.

## Generated outputs on `main`

Decided in #206, extended in #515. `main` carries code, tests, methods notes,
the small summary tables those notes quote numbers from, and the figure
outputs under `trait-visualization/`. Everything else is not a store for bulk
analysis output. A "note" here is a file under `methods/` or `docs/`.

- **Figures under `trait-visualization/`** are committed: PNG, SVG, HTML and
  PDF figures, the TSV and JSON (including graph JSON) they are drawn from,
  and Markdown write-ups. New files must be written by a script or `just`
  recipe in this repository, so a rerun replaces them in place. If a figure
  needs non-default script arguments, record the exact command in a README
  beside it. `archive/` subfolders hold frozen snapshots kept for comparison;
  nothing regenerates them. Don't commit shortcuts, caches or other local
  files. `trait_profile_investigation_figures.pptx` is a hand-built deck that
  a script only updates (#519).
- **Size limits** (`check-added-large-files` in pre-commit): 500 KB by
  default; 2 MiB for six-trait heatmap PNGs under
  `trait-visualization/trait-profile-investigation/six-shared-traits/` (#508).
  The default hook checks only files a commit adds, so it does not check
  merges or files already tracked; the six-trait PNG hook uses
  `--enforce-all` and also checks those. The figures that arrived when
  `clustered-heatmap-generation` landed (#515), 87 of them over 500 KB, were
  accepted as they are; new figures over the limit need a decision in the PR.
- **Commit** any other generated file only if a methods note cites it, it is
  under 500 KB, and it is regenerated by a `just` recipe named in that note.
- **Opt in path by path.** `data/derived/*` stays ignored; each committed file
  gets its own `!` line in `.gitignore`, so every addition shows in review. Do
  not `git add -f`. Git never looks inside an ignored directory, so a file in a
  subdirectory needs each level re-included and its other contents re-ignored:

  ```gitignore
  !data/derived/calibration/
  data/derived/calibration/*
  !data/derived/calibration/composition/
  data/derived/calibration/composition/*
  !data/derived/calibration/composition/composition_calibration_summary.tsv
  ```

  Check with `git check-ignore -v <path>` (no output means it is tracked).
- **Outside `trait-visualization/`, don't commit** raw null draws
  (per-replicate summaries are fine), caches, images, HTML, graph JSON or
  before/after snapshots. The note names the recipe that regenerates them. If a
  fixed copy has to be citable, attach it to a GitHub release (or Zenodo) and
  record its SHA-256 in the note.
- **No Git LFS.** `data/derived/study_trait_profiles_wide.tsv.gz` predates this
  rule and stays, checked against a 5 MiB limit (#508).

## Useful tools

### yq (recommended)

[yq](https://github.com/mikefarah/yq) (Mike Farah's Go version) is
useful for exploring the YAML output:

```bash
# macOS
brew install yq

# Linux
sudo wget -qO /usr/local/bin/yq https://github.com/mikefarah/yq/releases/latest/download/yq_linux_amd64
sudo chmod +x /usr/local/bin/yq

# Show all taxon names in a study
yq '.experiments[].signatures[].taxa[].taxon_name' data/studies/bsdb_10202341.yaml

# Count experiments
yq '.experiments | length' data/studies/bsdb_10202341.yaml

# Find studies with a specific PMID
grep -rl "pmid: 10202341" data/studies/
```

### jq (for JSON logs)

```bash
# Summarize ingest log skip reasons
jq '[.[] | select(.level == "skip")] | group_by(.reason) | map({reason: .[0].reason, count: length}) | sort_by(-.count)' data/derived/ingest_log.json
```

## Project conventions

- **File types**: CSV (input), YAML (output/schema), Python (code), JSON (logs).
  Nothing else.
- **Identifiers**: all entities have an `id` slot inherited from `NamedThing`.
  Study IDs use `bsdb:` prefix. Taxon IDs use `NCBITaxon:` prefix.
- **No sentinels**: missing data is `null`, not a magic value.
- **Log everything**: every data drop or transformation assumption goes
  to `data/derived/ingest_log.json`.
