#!/usr/bin/env bash
# merge_when_green.sh: merge a pull request only once every CI check passes
# (CLAUDE.md D52). Stands in for GitHub branch protection, which needs GitHub
# Pro on a private repository (the repo is public now, D72, and main is
# protected; this stays as a convenience: wait, then merge). Waits for the checks, stops on the first
# failure, and merges (merge commit, branch deleted) only if all passed.
#
# Usage:  bash .githooks/merge_when_green.sh <PR number>
set -euo pipefail
pr=${1:?Usage: merge_when_green.sh <PR number>}
repo=MacaulayWinter/seq-triage
echo "Waiting for CI on PR #$pr ..."
# GitHub registers a new commit's checks a few seconds after the push; asked
# too early, gh reports "no checks" and --watch gives up (it did twice,
# 2026-10-22). So wait, up to 5 minutes, until at least one check exists.
for _ in $(seq 1 30); do
    [ -n "$(gh pr checks "$pr" -R "$repo" 2>/dev/null)" ] && break
    sleep 10
done
# --watch waits until every check finishes; it exits non-zero if any failed.
if ! gh pr checks "$pr" -R "$repo" --watch --fail-fast --interval 30; then
    echo "CI did not pass on PR #$pr: not merged." >&2
    exit 1
fi
gh pr merge "$pr" -R "$repo" --merge --delete-branch
echo "PR #$pr merged."
