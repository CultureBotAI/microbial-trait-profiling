#!/usr/bin/env bash
# Purpose: show that the label-support co-assertion "both" count
#   (`either, both = ...` in src/mtp/label_support_audit.py) is caught by the
#   current tests (#1345, #1352), and which mutants the tests at 1bd7fa1 and
#   e9c0d36 missed. Only tracked files are copied, as they are in the working
#   tree (uncommitted edits to tracked files included; untracked and ignored files,
#   such as data directories and .venv, never), into a new directory under a
#   scratch folder outside the repository. Each mutant is applied there only, and
#   the label-support tests are run against the copy.
# Inputs: $1 = an existing scratch folder outside the repository; the script
#   refuses a folder inside it. Run from anywhere inside the repository: the root
#   comes from `git rev-parse --show-toplevel`. Uses the repository's .venv python.
#   Needs commits 1bd7fa1 and e9c0d36 in the history.
# Outputs: a new mutant-copy.* directory under $1, and one line per mutant and
#   test file on stdout: "<mutant> <test file> <pytest summary>".
# Issue: #1345 and #1352 (reviews of PR #1338 for #1237). Date: 2026-11-02.
set -euo pipefail

repo="$(git rev-parse --show-toplevel)"
repo_real="$(cd "$repo" && pwd -P)"
python="$repo/.venv/bin/python"
if [ "$#" -ne 1 ] || [ ! -d "$1" ]; then
	echo "usage: $0 <existing scratch folder outside the repository>" >&2
	exit 2
fi
scratch_real="$(cd "$1" && pwd -P)"
case "$scratch_real/" in
"$repo_real"/*)
	echo "refusing: $scratch_real is inside the repository $repo_real" >&2
	exit 2
	;;
esac

copy="$(mktemp -d "$scratch_real/mutant-copy.XXXXXX")"
git -C "$repo" ls-files -z | rsync -a --from0 --files-from=- "$repo/" "$copy/"
for commit in 1bd7fa1 e9c0d36; do
	git -C "$repo" show "$commit:tests/test_label_support_audit.py" \
		>"$copy/tests/test_label_support_audit_at_$commit.py"
done

target="src/mtp/label_support_audit.py"
original="either, both = a.gt(0) | b.gt(0), a.gt(0) & b.gt(0)"

apply_mutant() {
	"$python" - "$copy/$target" "$original" "$1" <<'EOF'
import sys
from pathlib import Path

path, old, new = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
text = path.read_text(encoding="utf-8")
if text.count(old) != 1:
    raise SystemExit(f"expected the original line exactly once in {path}")
path.write_text(text.replace(old, new), encoding="utf-8")
EOF
}

for spec in \
	"none|$original" \
	"both_never|either, both = a.gt(0) | b.gt(0), a.gt(0) & b.lt(0)" \
	"both_is_either|either, both = a.gt(0) | b.gt(0), a.gt(0) | b.gt(0)" \
	"both_needs_two|either, both = a.gt(0) | b.gt(0), a.gt(0) & b.gt(1)" \
	"both_is_side_a|either, both = a.gt(0) | b.gt(0), a.gt(0)" \
	"both_is_side_b|either, both = a.gt(0) | b.gt(0), b.gt(0)"; do
	name="${spec%%|*}"
	line="${spec#*|}"
	cp "$repo/$target" "$copy/$target"
	apply_mutant "$line"
	for tests in \
		tests/test_label_support_audit.py \
		tests/test_label_support_audit_at_1bd7fa1.py \
		tests/test_label_support_audit_at_e9c0d36.py; do
		summary="$(cd "$copy" && PYTHONPATH="$copy/src" "$python" -m pytest -q -p no:cacheprovider \
			"$tests" 2>&1 | tail -1 || true)"
		echo "$name $tests $summary"
	done
done
cp "$repo/$target" "$copy/$target"
