#!/usr/bin/env bash
# merge_when_green.sh: merge a pull request only once every CI check passes
# (CLAUDE.md D52). Stands in for GitHub branch protection, which needs GitHub
# Pro on a private repository. Waits for the checks, stops on the first
# failure, and merges (merge commit, branch deleted) only if all passed.
#
# Usage:  bash .githooks/merge_when_green.sh <PR number>
set -euo pipefail
pr=${1:?Usage: merge_when_green.sh <PR number>}
repo=MacaulayWinter/seq-triage
echo "Waiting for CI on PR #$pr ..."
# --watch waits until every check finishes; it exits non-zero if any failed.
if ! gh pr checks "$pr" -R "$repo" --watch --fail-fast --interval 30; then
    echo "CI did not pass on PR #$pr: not merged." >&2
    exit 1
fi
gh pr merge "$pr" -R "$repo" --merge --delete-branch
echo "PR #$pr merged."
