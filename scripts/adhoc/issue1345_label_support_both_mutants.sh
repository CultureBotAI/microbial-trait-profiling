#!/usr/bin/env bash
# Purpose: show that the label-support co-assertion "both" count
#   (`either, both = ...` in src/mtp/label_support_audit.py) is now caught by a
#   test (#1345), and that the tests at 1bd7fa1 did not catch it. The checkout is
#   copied to a new directory under a scratch folder outside the repository; each
#   mutant is applied there only, and the label-support tests are run against it.
# Inputs: $1 = scratch folder outside the repository. Run from the repository
#   root; uses the repository's .venv python. Needs commit 1bd7fa1 in the history.
# Outputs: a new mutant-copy.* directory under $1, and one line per mutant and
#   test file on stdout: "<mutant> <test file> <pytest summary>".
# Issue: #1345 (review of PR #1338 for #1237). Date: 2026-11-02.
set -euo pipefail

repo="$(pwd)"
python="$repo/.venv/bin/python"
scratch="$1"
case "$(cd "$scratch" && pwd)/" in
"$repo"/*) echo "the scratch folder must be outside the repository" >&2; exit 2 ;;
esac
copy="$(mktemp -d "$scratch/mutant-copy.XXXXXX")"
rsync -a --exclude .venv --exclude .git --exclude '__pycache__' "$repo/" "$copy/"
git -C "$repo" show 1bd7fa1:tests/test_label_support_audit.py \
	>"$copy/tests/test_label_support_audit_at_1bd7fa1.py"

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
	"both_needs_two|either, both = a.gt(0) | b.gt(0), a.gt(0) & b.gt(1)"; do
	name="${spec%%|*}"
	line="${spec#*|}"
	cp "$repo/$target" "$copy/$target"
	apply_mutant "$line"
	for tests in tests/test_label_support_audit.py tests/test_label_support_audit_at_1bd7fa1.py; do
		summary="$(cd "$copy" && PYTHONPATH="$copy/src" "$python" -m pytest -q -p no:cacheprovider \
			"$tests" 2>&1 | tail -1 || true)"
		echo "$name $tests $summary"
	done
done
cp "$repo/$target" "$copy/$target"
